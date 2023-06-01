#include "ElasticBody.h"

extern const uint8_t UNIFORM, LOG, DYNAMIC;
extern const uint8_t FLAT, POLY, SPHERE;
extern const geometry HERTZ;

double ellipe(double m) { 
/* 
  complete elliptic integral of the 2nd kind with parameter m = k^2 < 1 in R.
  The case m<0 is an implementation of eqs. (17.4.17) and (17.4.18) in the
  Handbook of Mathematical Functions by Abramowitz and Stegun.
  See also: https://math.stackexchange.com/questions/2461844/incomplete-elliptic-integral-of-the-second-kind-with-negative-parameter
*/

  using boost::math::ellint_1, boost::math::ellint_2, boost::math::jacobi_sn, boost::math::jacobi_cd;
  if (m>=0) return ellint_2(sqrt(m));

  double M = -m;
  double s = sqrt(1.+M), sqrt_b = sqrt(M/(1.+M));
  double u = ( ellint_1(sqrt_b) - ellint_1(sqrt_b, 0) )/s; // = K(-m), see (17.4.17)
  double a = u*s;
  double sn = jacobi_sn(sqrt_b, a), cd = jacobi_cd(sqrt_b, a);
  double phi = asin(sn); // Jacobi Amplitude, see https://dlmf.nist.gov/22.16
  return s*ellint_2(sqrt_b, phi) - (M/s)*sn*cd; // = E(-m), see (17.4.18)
}



double s_ana(double x) {
  if (abs(x-1) < 1.e-10) return 0;
  return ellipe(-4*x/pow(x-1, 2))/(M_PI*pow(x+1, 2)*abs(x-1));
};



ElasticBody::ElasticBody(vector<double>& edges, vector<double>& centers, 
                         map<string,string> elastic) : ID(0) {
  init(edges, centers, elastic);
}

ElasticBody::ElasticBody(vector<double>& edges, vector<double>& centers, 
                         map<string,string> elastic, uint32_t i) : ID(i) {
  init(edges, centers, elastic);
}



void ElasticBody::init(vector<double>& edges, vector<double>& centers, map<string,string> toml) {
  cerr << "# ElasticBody::init()\n";//FLOW

  initialized = 0;

  for (const auto& [key, value] : toml) {
    //cout << key << " " << value << endl;//DEBUG
    size_t idx = key.find(".");
    if (idx != string::npos && key.substr(0,idx) != "ElasticBody") continue;//ADJUST for ID
    string name = key.substr(idx+1);
    if (name == "Estar") Estar = stod(value);
    // solver
    else if (name == "damping") damping = stod(value);
    else if (name == "massGFMD") massGFMD = stod(value);
    // viscoelasticity
    else if (name == "nMaxwell") nMaxwell = stoi(value);
    else if (name == "Estar_inf") Estar_inf = stod(value);
    // force/displacement control
    else if (name == "force") force = stod(value);
    else if (name == "fConstraint") fConstraint = stoi(value);
    else if (name == "zPos") zPos = stod(value);
    //else std::cout << "# unknown elastic parameter: " << name << std::endl;//DEBUG
  }

  nBin = centers.size();
  bin_edge.resize(nBin+1);
  bin_center.resize(nBin);
  area = M_PI*edges.back()*edges.back();
  for (int iEdge = 0; iEdge <= nBin; ++iEdge) 
    bin_edge[iEdge] = edges[iEdge];
  for (int iBin = 0; iBin < nBin; ++iBin) 
    bin_center[iBin] = 0.5*(bin_edge[iBin+1] + bin_edge[iBin]);
  initialized += 1;  
  
  init_fields_from_bins();
  
  if (nMaxwell > 0) init_Maxwell(stod(toml["dTime"]));
  else Estar_inf = Estar;

  init_stiffness();
  moni.open("moni"+to_string(ID)+".dat");
  //moni = fopen("moni.dat", "w+");
}



void ElasticBody::init_fields_from_bins() {
/* 
  sets up all other fields according to the generated bin edges
*/
  uint32_t nBin = bin_edge.size() - 1;

  bin_area.resize(nBin);
  ext_stress.resize(nBin);
  int_stress.resize(nBin);
  disp.resize(nBin);
  disp_old.resize(nBin);
  stiffness_array.resize(nBin);

  for (int iBin = 0; iBin < nBin; ++iBin) {
    stiffness_array[iBin].resize(nBin);
    bin_area[iBin] = 2*M_PI*bin_center[iBin]*(bin_edge[iBin+1] - bin_edge[iBin]);
  }

  initialized += 2;
};



void ElasticBody::init_stiffness() {
/* 
  stiffness represents distribution of stress per surface displacement when 
  only a single circle of width <bin_width> at r=<bin_center> is displaced.
*/
  cerr << "# ElasticBody::init_stiffness()\n";//FLOW

  double max_stiff_loc = -1e256;
  for (int uBin = 0; uBin < nBin; ++uBin) {
    double bin_width = bin_area[uBin]/(2*M_PI*bin_center[uBin]);
    double prefac = -Estar_inf*bin_width/pow(bin_center[uBin], 2);
    for (int sBin = 0; sBin < nBin; ++sBin) {
      if (uBin==sBin) continue;
      double x = bin_center[sBin]/bin_center[uBin];
      stiffness_array[sBin][uBin] = prefac*s_ana(x);
    }
  }

  // center of mass displacement must not result in any stress
  for (int sBin = 0; sBin < nBin; ++sBin) {
    double stress_per_u = 0;
    for (int uBin = 0; uBin < nBin; ++uBin) {
      if (uBin==sBin) continue;
      stress_per_u -= stiffness_array[sBin][uBin];
    }
    stiffness_array[sBin][sBin] = stress_per_u;
    if (stiffness_array[sBin][sBin] > max_stiff_loc) max_stiff_loc = stiffness_array[sBin][sBin];
  }
  //cout << "max stiffness: " << max_stiff_loc << "\n";//DEBUG*/
  max_stiff = (Estar/Estar_inf)*max_stiff_loc;

  initialized += 4;
};



void ElasticBody::init_Maxwell(double dTime){
  cerr << "# ElasticBody::init_Maxwell()\n";//FLOW
  dispMw.resize(nMaxwell);
  invTauMw.resize(nMaxwell);
  stiffFacMw.resize(nMaxwell);
  for (vector<vector<double> >& array : dispMw) {
    array.resize(bin_center.size());
    for (vector<double>& vec : array) {
      vec.resize(bin_center.size());
      for (double& val : vec) val = 0;
    }
  }

  string filename = "maxwell.in";
  ifstream input(filename);
  int nMaxwell_test;
  double dTime_test=0, kLow_test=0, kHigh_test=0, damping_test=0, mass_test=0;

  if(!input.is_open()) terminate(filename+" does not exist.");

  cerr << "# Reading "+filename+".\n";

  size_t pos; // current position in file stream
  std::string ROL; //REST OF LINE
  std::size_t NIS = std::string::npos; //NIS = NOT IN STRING
  
  while (input.peek()!=EOF){

    // skip empty lines and lines starting with '#'
    char firstChar;
    pos = input.tellg();
    input >> firstChar;
    input.seekg(pos,input.beg);
    if (firstChar=='#' || firstChar=='\n') {
      getline(input, ROL);
      continue;
    }

    // read stiffness, stiffHigh, nMaxwell
    pos = input.tellg();
    double param;
    input >> param;
    getline(input, ROL);
    if (ROL.find("# stiffness #") !=NIS) {
      kLow_test = param;
      if (kLow_test != Estar) terminate("stiffness0 incompatible with "+filename);
    }
    if (ROL.find("# stiffHigh #") !=NIS) {
      kHigh_test = param;
      if (kHigh_test != Estar_inf) terminate("stiffHigh0 stiffness incompatible with "+filename);
    }
    if (ROL.find("# nMaxwell #") !=NIS) {
      nMaxwell_test = param;
      if (nMaxwell_test != nMaxwell) terminate("nMaxwell incompatible with "+filename);
    }
    if (ROL.find("# dTime #") !=NIS) dTime_test = param;
    if (ROL.find("# massGFMD #") !=NIS) {
      mass_test = param;
      if (mass_test != massGFMD) terminate("massGFMD incompatible with "+filename);
    }
    if (ROL.find("# dampGlobal #") !=NIS) damping_test = param;
    if (ROL.find("# MWelement #") !=NIS) {
      double iMw=-1, stiffMw=-1, tauMw=-1;

      input.seekg(pos,input.beg);
      input >> iMw >> stiffMw >> tauMw;
      getline(input,ROL);
      invTauMw[iMw] = 1./tauMw;
      stiffFacMw[iMw] = stiffMw/Estar_inf;
    }
  }

  // sanity checks
  if (dTime_test < 0.99*dTime) {
    cerr << "CAUTION! dTime likely too large. Consider using this one or smaller:\n";
    cerr << dTime_test << "\t\t# dTime #\n";
  } 
  if (damping_test != damping) cerr << "CAUTION! You are using a different dampGlobal than suggested by maxwell.cpp!\n";

  initialized += 8;
  //cerr << "Finished reading maxwell.in.\n";//DEBUG
}



void ElasticBody::set_stress(double val) {
  for (int i = 0; i < disp.size(); ++i) {
    int_stress[i] = val;
    ext_stress[i] = val;
  }
};

void ElasticBody::set_disp(double val) {
  for (int i = 0; i < disp.size(); ++i) disp[i] = val;
};

void ElasticBody::set_values(vector<double> source, vector<double> *target) {
  if (source.size() != (*target).size()) 
    terminate("interpolation not implemented yet in ElasticBody::set_values().");
  for (int i = 0; i < source.size(); ++i) {
    (*target)[i] = source[i];
  }
};

void ElasticBody::set_disp_old(double val) {
  for (int i = 0; i < disp.size(); ++i) disp_old[i] = val;
};

void ElasticBody::set_damping(double val) {damping = val;};



void ElasticBody::stress_internal() {
  for (int sBin = 0; sBin < int_stress.size(); ++sBin) {
    int_stress[sBin] = 0;
    for (int uBin = 0; uBin < disp.size(); ++uBin) {
      int_stress[sBin] -= disp[uBin]*stiffness_array[sBin][uBin];
    }
  }

  /*// add Maxwell stress (local and non-local viscoelasticity)
  for (int iMw = 0; iMw < nMaxwell; ++iMw) {
    for (uint32_t sBin = 0; sBin < bin_center.size(); ++sBin) {
      double stressMw = 0;
      for (int iMw = 0; iMw < nMaxwell; ++iMw) {
        for (uint32_t uBin = 0; uBin < bin_center.size(); ++uBin) {
          stressMw += stiffFacMw[iMw]*stiffness_array[sBin][uBin]*dispMw[iMw][sBin][uBin];
          //stressMw += stiffFacMw[iMw]*abs(stiffness_array[sBin][uBin])*dispMw[iMw][sBin][uBin];
        }
      }
      int_stress[sBin] += stressMw;
    }
  }//*/

  // add Maxwell stress (purely local viscoelasticity)
  for (int iMw = 0; iMw < nMaxwell; ++iMw) {
    for (uint32_t sBin = 0; sBin < bin_center.size(); ++sBin) {
      int_stress[sBin] += stiffFacMw[iMw]*stiffness_array[sBin][sBin]*dispMw[iMw][sBin][sBin];
    }
  }//*/
};



void ElasticBody::stress_external() {
  double pressure = force/area;
  for (int iBin = 0; iBin < nBin; ++iBin) {
    ext_stress[iBin] += 1.*pressure;
  }
};



void ElasticBody::propagate(double dTime) {
  //cout << "dTime=" << dTime << "\n";//DEBUG
  //cout << "damping=" << damping << "\n";//DEBUG
  double dt2 = dTime*dTime;
  double damp_eff = 0.5*damping*dTime;

  double weighted_sum = 0;
  for (int iBin = 0; iBin < disp.size(); ++iBin) {
    // Langtangen Verlet (http://hplgit.github.io/fdm-book/doc/pub/vib/html/._vib003.html)
    //double inv_mass = max_stiff;//TEMP
    double inv_mass = stiffness_array[iBin][iBin]/massGFMD;
    double disp_new = 2*disp[iBin] + (damp_eff-1)*disp_old[iBin];
    //disp_new += (int_stress[iBin]+ext_stress[iBin])*dt2*inv_mass;//stress
    disp_new += (int_stress[iBin]+ext_stress[iBin])*bin_area[iBin]*dt2*inv_mass;//force
    disp_new /= 1 + damp_eff;

    // step forward
    disp_old[iBin] = disp[iBin];
    disp[iBin] = disp_new;
    weighted_sum += disp_new*bin_area[iBin];
  }
  dispCOM = weighted_sum/area;

  // update Maxwell elements via improved Euler's method (a.k.a. Heun's method)
  for (int iMw = 0; iMw < nMaxwell; ++iMw) {
    for (uint32_t uBin = 1; uBin < bin_center.size(); ++uBin) {
      for (uint32_t sBin = 1; sBin < bin_center.size(); ++sBin) {
        double slopeNow = invTauMw[iMw]*(disp[uBin]-dispCOM - dispMw[iMw][sBin][uBin]);
        double dispFnew = dispMw[iMw][sBin][uBin] + slopeNow*dTime;
        double slopeNew = invTauMw[iMw]*(disp[uBin]-dispCOM - dispFnew);
        dispMw[iMw][sBin][uBin] += 0.5*(slopeNow + slopeNew) * dTime;
      }
    }
  }

  moni << disp[0] << "\t" << disp.back() << "\n";
  //fprintf(moni, "%f\t%f\n", disp[0], disp.back());
  //printf("%f\t%f\n", disp[0], disp.back());
};

void ElasticBody::test_ellip_int() {
  cerr << "# ElasticBody::test_ellip_int()\n";//FLOW
  /*
  outputs value pairs (m, ellipe(m)). This way, the present implementation 
  can directly be compared against scipy.special's ellipe(m).
  */

  double x = -1.99, dx = 0.01;
  ofstream output("elliptic2.dat");
  output << "# k^2\tE(k^2)\n";
  while (x < 1.) {
    output << x;
    output << "\t" << ellipe(x);
    output << "\n";
    x += dx;
  }
  output.close();
}

void ElasticBody::test_Verlet() {
  cerr << "# test_Verlet()\n";//FLOW
  bin_edge = {0,2};
  bin_center = {1};
  nBin = 1;
  init_fields_from_bins();
  init_stiffness();
  stiffness_array[0][0] = 1.;
  set_disp(1.); // u(t=0) relaxation
  set_disp_old(0.999);
  set_damping(0.0);

  ofstream output("verlet.dat");
  output << "# time\tdisplacement\n";
  double dTime = 0.1;
  uint32_t nTime = 1000;
  for (int iTime = 0; iTime < nTime; ++iTime) {
    set_stress(0);
    stress_internal();
    propagate(dTime);
    output << iTime*dTime << "\t" << disp[0] << "\n";
  }
  output.close();
}



void ElasticBody::write_params(const string& filename) {
  ofstream output(filename, std::ios::app);
  if (!output.is_open()) return;

  output << "[ElasticBody]\n";
  output << "  force = " << force << "\n";
  output << "  Estar = " << Estar << "\n";
  output << "  nMaxwell = " << nMaxwell << "\n";
  if (nMaxwell) output << "  Estar_inf = " << Estar_inf << "\n";
  output << "  damping = " << damping << "\n";
  output << "  massGFMD = " << massGFMD << "\n";
  output << "# max_stiff = " << get_stiff() << "\n";
  output << "\n";

  output.close();
}


void ElasticBody::write_config() {
  io::write_array("stiffness.dat", stiffness_array, "v sBin v || < uBin >");
  io::write_vectors("config.dat", 
      {&bin_center, &disp, &int_stress, &ext_stress}, "r\tdisp\tint_stress\text_stress");
}

void ElasticBody::write_config(uint32_t iTime) {
  io::write_vectors("config."+to_string(iTime)+".dat", 
      {&bin_center, &disp, &int_stress, &ext_stress}, "r\tdisp\tint_stress\text_stress");
}
