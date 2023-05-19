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


void ElasticBody::init(uint32_t N, double Rmax, uint8_t kind) {
  if (kind==UNIFORM) init_uniform_bins(N, Rmax);
  else {
    cerr << "non-uniform grids not implemented yet.\n";
    exit(1);
  }
  init_stiffness();
}

ElasticBody::ElasticBody(uint32_t N, double Rmax, uint8_t kind) {
  init(N, Rmax, kind);
}

ElasticBody::ElasticBody(map<string,string> global, map<string,string> elastic) {
  uint32_t nBin;
  double Rmax;
  uint8_t grid;
  for (const auto& [key, value] : global) {
    cout << key << " " << value << endl;//DEBUG
    if (key == "Rmax") Rmax = stod(value);
    else if (key == "nBin") nBin = stoi(value);
    else if (key == "grid") grid = stoi(value);
    else std::cout << "# unknown global parameter: " << key << std::endl;//DEBUG
  }

  for (const auto& [key, value] : elastic) {
    cout << key << " " << value << endl;//DEBUG
    if (key == "Estar") Estar = stod(value);
    // solver
    else if (key == "damping") damping = stod(value);
    else if (key == "massGFMD") massGFMD = stod(value);
    // viscoelasticity
    else if (key == "nMaxwell") nMaxwell = stoi(value);
    else if (key == "Estar_inf") Estar_inf = stod(value);
    // force/displacement control
    else if (key == "force") force = stod(value);
    else if (key == "fConstraint") fConstraint = stoi(value);
    else if (key == "zPos") zPos = stod(value);
    else std::cout << "# unknown elastic parameter: " << key << std::endl;//DEBUG
  }

  if (nMaxwell) {
    dispMw.resize(nMaxwell);
    for (vector<double> v : dispMw) v.resize(nBin);
  }
  else Estar_inf = Estar;

  init(nBin, Rmax, grid);
  if (nMaxwell) init_Maxwell(stod(global["dTime"]));//TODO: improve

}

void ElasticBody::init_uniform_bins(uint32_t N, double Rmax) {
/* 
  discretizes radial coordinate into <N> equidistant points from 0 to <Rmax>. 
*/

  nBin = N;
  bin_edge.resize(nBin+1);
  area = M_PI*Rmax*Rmax;

  double dr = Rmax/nBin;
  for (int iEdge = 0; iEdge <= nBin; ++iEdge) 
    bin_edge[iEdge] = iEdge*dr;
  
  init_fields_from_bins();
};

void ElasticBody::init_fields_from_bins() {
/* 
  sets up all other fields according to the generated bin edges
*/

  bin_center.resize(nBin);
  bin_area.resize(nBin);
  ext_stress.resize(nBin);
  int_stress.resize(nBin);
  disp.resize(nBin);
  disp_old.resize(nBin);
  stiffness_array.resize(nBin);

  for (int iBin = 0; iBin < nBin; ++iBin) {
    stiffness_array[iBin].resize(nBin);
    bin_center[iBin] = 0.5*(bin_edge[iBin+1] + bin_edge[iBin]);
    bin_area[iBin] = 2*M_PI*bin_center[iBin]*(bin_edge[iBin+1] - bin_edge[iBin]);
  }

  initialized = 1;
};


void ElasticBody::init_stiffness() {
/* 
  stiffness represents distribution of stress per surface displacement when 
  only a single circle of width <bin_width> at r=<bin_center> is displaced.
*/
  //cerr << "# ElasticBody::init_stiffness()\n";//FLOW

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
  max_stiff = max_stiff_loc;

  initialized = 2;
};

void ElasticBody::init_Maxwell(double dTime){
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

  //cerr << "Finished reading maxwell.in.\n";//DEBUG
}


void ElasticBody::read_stiffness(const string& filename) {
  //cerr << "# ElasticBody::read_stiffness()\n";//FLOW
  ifstream input(filename);
  if (!input.is_open()) {
    cerr << filename+" not found.\n";
    exit(1);
  }
  cerr << "# reading "+filename+"\n";
  string skip_string; double skip_value;
  for (int sBin = 0; sBin < nBin; ++sBin) {
    if (input.peek() == '#') getline(input, skip_string);
    //input >> skip_value; // radius
    for (int uBin = 0; uBin < nBin; ++uBin) input >> stiffness_array[sBin][uBin];
    getline(input, skip_string);
  }
  input.close();
};

void ElasticBody::set_stress(double val) {
  for (int i = 0; i < disp.size(); ++i) {
    int_stress[i] = val;
    ext_stress[i] = val;
  }
};

void ElasticBody::set_disp(double val) {
  for (int i = 0; i < disp.size(); ++i) disp[i] = val;
};

void set_values(vector<double> source, vector<double> *target) {
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

void ElasticBody::internal_stress() {
  for (int sBin = 0; sBin < int_stress.size(); ++sBin) {
    int_stress[sBin] = 0;
    for (int uBin = 0; uBin < disp.size(); ++uBin) {
      int_stress[sBin] -= disp[uBin]*stiffness_array[sBin][uBin];
    }
  }

  for (int iMw = 0; iMw < nMaxwell; ++iMw) {
    terminate("Maxwell not properly implemented yet.");

    // add Maxwell stress
    for (uint32_t iBin = 0; iBin < nBin; ++iBin) {
      double stressMw = 0;
      for (int iMw = 0; iMw < nMaxwell; ++iMw) {
        stressMw += stiffFacMw[iMw]*stiffness_array[iBin][iBin]*dispMw[iMw][iBin];
      }
      int_stress[iBin] += stressMw;//TEMP isolate conventional GFMD element by commenting this out
    }
  }
};

void ElasticBody::external_stress() {
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

  for (int iBin = 0; iBin < disp.size(); ++iBin) {
    // Langtangen Verlet (http://hplgit.github.io/fdm-book/doc/pub/vib/html/._vib003.html)
    double inv_mass = stiffness_array[iBin][iBin]/massGFMD;
    double disp_new = 2*disp[iBin] + (damp_eff-1)*disp_old[iBin];
    //disp_new += (int_stress[iBin]+ext_stress[iBin])*dt2*inv_mass;
    disp_new += (int_stress[iBin]+ext_stress[iBin])*bin_area[iBin]*dt2*inv_mass;
    disp_new /= 1 + damp_eff;

    // step forward
    disp_old[iBin] = disp[iBin];
    disp[iBin] = disp_new;
  }

  for (int iMw = 0; iMw < nMaxwell; ++iMw) {
    terminate("Maxwell not properly implemented yet.");
    // update Maxwell elements via improved Euler's method (a.k.a. Heun's method)
    for (uint32_t iBin = 1; iBin < nBin; ++iBin) {
      for (int iMw = 0; iMw < nMaxwell; ++iMw) {
        double slopeNow = invTauMw[iMw]*(disp[iBin] - dispMw[iMw][iBin]);
        double dispFnew = dispMw[iMw][iBin] + slopeNow*dTime;
        double slopeNew = invTauMw[iMw]*(disp[iBin] - dispFnew);
        dispMw[iMw][iBin] += 0.5*(slopeNow + slopeNew) * dTime;
      }
    }
  }
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





// -------------------------------------------------------------------------- //



double s_fit(double x) {
  return 1./(pow(x,2)-pow(abs(x),3)/M_PI);
};

void ElasticBody::fit_stiffness() {
/* 
  stiffness represents distribution of stress per surface displacement when 
  only a single circle of width <bin_width> at r=<bin_center> is displaced.
  This is an old version of init_stiffness, where eqs were fitted to GFMD data.
*/
  //cerr << "# ElasticBody::fit_stiffness()\n";//FLOW
  
  double max_stiff = -1e256;
  for (int uBin = 0; uBin < nBin; ++uBin) {

    double bin_width = bin_area[uBin]/(2*M_PI*bin_center[uBin]);
    double prefac = -Estar*bin_width/(2*M_PI*bin_center[uBin]*bin_center[uBin]);

    double Fout_per_u = 0; // total force acting outside the displaced circle
    for (int sBin = 0; sBin < nBin; ++sBin) {
      if (uBin==sBin) continue;

      double x = bin_center[sBin]/bin_center[uBin];
      
      /*// OPTION: separate treatment for displacement at the center
      //if (uBin==0) stiffness_array[sBin][uBin] = M_PI*prefac/(pow(x-1., 3)/M_PI);
      if (uBin==0) stiffness_array[sBin][uBin] = 2.4*prefac/(pow(x-1., 3)/M_PI);//works surprisingly well*/
      
      // OPTION: combination of the two below
      if (x>1) stiffness_array[sBin][uBin] = prefac/(pow(x-1., 2) + pow(x-1., 3)/M_PI);
      //else stiffness_array[sBin][uBin] = prefac*(s_fit(x-1)+s_fit(x+1) + 0.2);
      //else stiffness_array[sBin][uBin] = prefac*(pow(x-1., -2) + pow(x+1., -2));//simple without constant
      else stiffness_array[sBin][uBin] = prefac*(pow(x-1., -2) + 3.0);//works surprisingly well*/

      /*// OPTION: x**-2 only
      if (x>1) stiffness_array[sBin][uBin] = prefac*pow(x-1., -2);
      else stiffness_array[sBin][uBin] = prefac*(pow(x-1., -2) + pow(x+1., -2) + M_PI/2.);//*/
      
      /*// OPTION: x**-2 to x**-3 cross-over
      if (x>1) stiffness_array[sBin][uBin] = prefac/(pow(x-1., 2) + pow(x-1., 3)/M_PI);
      else stiffness_array[sBin][uBin] = prefac*(1./(pow(x-1., 2) + pow(x+1., 2)) + 2*M_PI/5.5);//*/
      
      Fout_per_u += stiffness_array[sBin][uBin]*bin_area[uBin];
    }

    // force inside the circle must balance out the force outside of it
    //stiffness_array[uBin][uBin] = -Fout_per_u/bin_area[uBin];
  }

  // // find constant term so that two sum rules can be fulfilled at the same time
  // for (int iBin = 0; iBin < nBin; ++iBin) {
  //   double norm = 0, offset = 0;
  //   for (int jBin = 0; jBin < nBin; ++jBin) {
  //     if (iBin==jBin) continue;
  //     double ratio = bin_area[jBin]/bin_area[iBin];
  //     offset += ratio*stiffness_array[jBin][iBin];
  //     offset -= stiffness_array[iBin][jBin];
  //     if (jBin>iBin) norm ++;
  //     else norm -= ratio;
  //   }
  //   double c_stiff = offset/norm;
  //   cout << c_stiff << "\n";//DEBUG
  // }

  // center of mass displacement must not result in any stress
  for (int sBin = 0; sBin < nBin; ++sBin) {
    double stress_per_u = 0;
    for (int uBin = 0; uBin < nBin; ++uBin) {
      if (uBin==sBin) continue;
      stress_per_u -= stiffness_array[sBin][uBin];
    }
    stiffness_array[sBin][sBin] = stress_per_u;
    //stiffness_array[sBin][sBin] = 1.00*stress_per_u - 0.00*stiffness_array[sBin][sBin];
    if (stiffness_array[sBin][sBin] > max_stiff) max_stiff = stiffness_array[sBin][sBin];
  }
  //cout << "max stiffness: " << max_stiff << "\n";//DEBUG*/
};