#include "header.h"
//#include "elast.cpp"
//#include "indenter.cpp"

// All DOFs are set to a circular eigenfrequency of sqrt(k/m) = 1.
// Therefore, 2.0 is exactly critical damping
double damping = 1.8;
double dTime = 0.05;

namespace io {
  template <typename T> vector<T> read_column(const string& filename, uint32_t col);
  double column_increment(const string&, uint32_t);
  void write_vectors(const string& filename, vector< vector<double>* > arrays, const string& header);
  void write_array(const string& filename, vector< vector<double> >& array, const string& header);
}

namespace elast {
  static double Estar = 1.;
  static double force = 0.1, area = M_PI*4;
  static uint32_t nBin = 32;
  static vector<double> bin_center(nBin), bin_area(nBin); // length nBin
  static vector<double> ext_stress(nBin), int_stress(nBin), disp(nBin), disp_old(nBin); // length nBin
  static vector<double> bin_edge; // length nBin+1
  static vector<vector<double> > stiffness_array;

  void init_uniform_bins(uint32_t N, double Rmax);
  void init_fields_from_bins();
  void init_stiffness();
  void fit_stiffness();
  void read_stiffness(const string&);

  void set_stress(double val);
  void set_disp(double val);
  void set_disp_old(double val);
  void internal_stress();
  void external_stress();
  void compute_stress();
  void propagate();

  void test_stress();
  void test_disp_propagation();
}//*/

namespace indenter {
  static double z0 = 0;
  static double dzInd = 1e-3;
  static double radius = 1.;
  static double exponent = 2.;
  static vector<double> height;

  static double polynomial(double r) {
    double h = pow(r*r/radius, exponent/2)/exponent;
    return(z0 + h);
  };

  static void init_height(vector<double>& rval, double (*zpos)(double)) {
    height.resize(rval.size());
    for (int ir = 0; ir < rval.size(); ++ir) {
      height[ir] = zpos(rval[ir]);
    }
  };
}//*/

namespace interaction {
  double curvature = 10.;
  double surfEnerg = 1e-4;
  double range = M_PI*sqrt(surfEnerg/(2*curvature));

  double potential1(double gap) {
    if (gap > 0) return 0;
    return curvature*gap;
  };

  double potential2(double gap) {
    if (gap > range) return 0;
    else if (gap > 0) return surfEnerg*M_PI/range*sin(gap*M_PI/range)/2.;
    return curvature*gap;
  };
}//*/



void elast::init_uniform_bins(uint32_t N, double Rmax) {
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

void elast::init_fields_from_bins() {
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
};


double ellipe(double m) { 
  /* 
  complete elliptic integral of the 2nd kind with parameter m = k^2 < 1 in R.
  The case m<0 is an implementation of eqs. (17.4.17) and (17.4.18) in the
  Handbook of Mathematical Functions by Abramowitz and Stegun.
  See also: https://math.stackexchange.com/questions/2461844/incomplete-elliptic-integral-of-the-second-kind-with-negative-parameter
  */

  using namespace boost::math;
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
  return -ellipe(-4*x/pow(x-1, 2))/(M_PI*pow(x+1, 2)*abs(x-1));
};

void elast::init_stiffness() {
  cerr << "# init_stiffness()\n";//FLOW
  /* 
  stiffness represents distribution of stress per surface displacement when 
  only a single circle of width <bin_width> at r=<bin_center> is displaced.
  */

  double max_stiff = -1e256;
  for (int uBin = 0; uBin < nBin; ++uBin) {
    double bin_width = bin_area[uBin]/(2*M_PI*bin_center[uBin]);
    double prefac = -Estar*bin_width/(2*M_PI*bin_center[uBin]*bin_center[uBin]);
    for (int sBin = 0; sBin < nBin; ++sBin) {
      if (uBin==sBin) continue;
      double x = bin_center[sBin]/bin_center[uBin];
      stiffness_array[sBin][uBin] = s_ana(x)*bin_width/pow(bin_center[uBin], 2);
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
    if (stiffness_array[sBin][sBin] > max_stiff) max_stiff = stiffness_array[sBin][sBin];
  }
  cout << "max stiffness: " << max_stiff << "\n";//DEBUG*/
};


double s_fit(double x) {
  return 1./(pow(x,2)-pow(abs(x),3)/M_PI);
};

void elast::fit_stiffness() {
  cerr << "# fit_stiffness()\n";//FLOW
  /* 
  stiffness represents distribution of stress per surface displacement when 
  only a single circle of width <bin_width> at r=<bin_center> is displaced.
  This is an old version of init_stiffness, where eqs were fitted to GFMD data.
  */

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
  cout << "max stiffness: " << max_stiff << "\n";//DEBUG*/
};

void elast::read_stiffness(const string& filename) {
  cerr << "# read_stiffness()\n";//FLOW
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

void elast::set_stress(double val) {
  for (int i = 0; i < disp.size(); ++i) {
    int_stress[i] = val;
    ext_stress[i] = val;
  }
};

void elast::set_disp(double val) {
  for (int i = 0; i < disp.size(); ++i) disp[i] = val;
};

void elast::set_disp_old(double val) {
  for (int i = 0; i < disp.size(); ++i) disp_old[i] = val;
};

void elast::internal_stress() {
  for (int sBin = 0; sBin < int_stress.size(); ++sBin) {
    int_stress[sBin] = 0;
    for (int uBin = 0; uBin < disp.size(); ++uBin) {
      int_stress[sBin] -= disp[uBin]*stiffness_array[sBin][uBin];
    }
  }
};

void elast::external_stress() {
  double pressure = force/area;
  for (int iBin = 0; iBin < nBin; ++iBin) {
    ext_stress[iBin] = interaction::potential1(indenter::height[iBin]-disp[iBin]);
    ext_stress[iBin] += 1.*pressure;
  }
};

void elast::compute_stress() {
  set_stress(0);
  external_stress();
  internal_stress();
}

void elast::propagate() {
  double dt2 = dTime*dTime;
  double damp_eff = 0.5*damping*dTime;

  for (int iBin = 0; iBin < disp.size(); ++iBin) {
    // Langtangen Verlet (http://hplgit.github.io/fdm-book/doc/pub/vib/html/._vib003.html)
    double inv_mass = stiffness_array[iBin][iBin];
    double disp_new = 2*disp[iBin] + (damp_eff-1)*disp_old[iBin];
    //disp_new += (int_stress[iBin]+ext_stress[iBin])*dt2*inv_mass;
    disp_new += (int_stress[iBin]+ext_stress[iBin])*bin_area[iBin]*dt2*inv_mass;
    disp_new /= 1 + damp_eff;

    // step forward
    disp_old[iBin] = disp[iBin];
    disp[iBin] = disp_new;
  }

  /*//HACK center of mass motion
  double ext_force = 0;
  for (int iBin = 0; iBin < disp.size(); ++iBin) {
    ext_force += ext_stress[iBin]*bin_area[iBin];
  }
  double delta_z = (force+ext_force)/(1.33*Estar*bin_edge.back()*M_PI);
  for (int iBin = 0; iBin < disp.size(); ++iBin) {
    disp[iBin] += delta_z;
  }//*/
};

void test_ellip_int() {
  cerr << "# test_ellip_int()\n";//FLOW
  /*
  outputs value pairs (m, ellipe(m)). This way, the present implementation 
  can directly be compared against scipy.special's ellipe(m).
  */

  using namespace boost::math;
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

void elast::test_stress() {
  cerr << "# test_stress()\n";//FLOW
  init_uniform_bins(200, 8.);
  init_stiffness();
  //read_stiffness("stiff.dat");
  vector<double> solution(nBin);

  // flat punch: calculate stress from displacement profile
  double a = 1., dist = force/(2.*Estar*a);
  for (int iBin = 0; iBin < nBin ; ++iBin) {
    if (bin_center[iBin] < a) {
      disp[iBin] = 0;
      solution[iBin] = (Estar*dist/M_PI)*pow(a*a - pow(bin_center[iBin],2), -0.5);
    }
    else {
      disp[iBin] = dist - (2*dist/M_PI)*asin(a/bin_center[iBin]);
      solution[iBin] = 0.;
    }
  }
  internal_stress();
  for (double& s : int_stress) s += force/(M_PI*pow(bin_edge.back(),2));

  io::write_vectors("stress_flat.dat", {&bin_center, &int_stress, &solution, &disp}, 
                    "r\tstress(sim)\tstress(ana)\tdisp(sim)");
  
  // Hertz: calculate stress from displacement profile
  double R = 1., ac = pow(0.75*R*force/Estar, 1./3);
  dist = ac*ac/R;
  for (int iBin = 0; iBin < nBin ; ++iBin) {
    if (bin_center[iBin] < ac) {
      disp[iBin] = pow(bin_center[iBin], 2)/(2*R);
      solution[iBin] = (2.*Estar/M_PI/R)*sqrt(ac*ac - pow(bin_center[iBin], 2));
    }
    else {
      disp[iBin] = dist - (ac*ac/M_PI/R)*((2-pow(bin_center[iBin]/ac,2))*asin(ac/bin_center[iBin]) + sqrt(pow(bin_center[iBin],2) - ac*ac)/ac);
      solution[iBin] = 0.;
    }
  }
  internal_stress();
  for (double& s : int_stress) s += force/(M_PI*pow(bin_edge.back(),2));

  io::write_vectors("config.dat", 
                        {&bin_center, &disp, &int_stress}, "r\tdisp\tint_stress");
  io::write_vectors("stress_hertz.dat", {&bin_center, &int_stress, &solution, &disp}, 
                    "r\tstress(sim)\tstress(ana)\tdisp(sim)");
};


void elast::test_disp_propagation() {
  cerr << "# test_disp_propagation()\n";//FLOW
  init_uniform_bins(200, 8.);
  init_stiffness();
  io::write_array("stiff.dat", elast::stiffness_array, "v sBin v || < uBin >");
  //read_stiffness("stiff.dat");
  //io::write_array("stiff_read.dat", stiffness_array, "v sBin v || < uBin >");//DEBUG
  indenter::init_height(bin_center, &(indenter::polynomial));
  vector<double> solution(nBin);

  // Hertz: calculate displacement from stress profile
  double R = 1., ac = pow(0.75*R*force/Estar, 1./3);
  double dist = ac*ac/R, pressure = force/area;
  for (int iBin = 0; iBin < nBin ; ++iBin) {
    if (bin_center[iBin] < ac) {
      ext_stress[iBin] = 1.*pressure - (2./M_PI/R)*sqrt(ac*ac - pow(bin_center[iBin], 2));
      solution[iBin] = pow(bin_center[iBin], 2)/(2*R);
    }
    else {
      ext_stress[iBin] = 1.*pressure;
      solution[iBin] = dist - (ac*ac/M_PI/R)*((2-pow(bin_center[iBin], 2)/(ac*ac))*asin(ac/bin_center[iBin]) + sqrt(pow(bin_center[iBin],2) -ac*ac)/ac);
    }
  }
  uint32_t nTime = 20000;
  for (int iTime = 0; iTime < nTime; ++iTime) {
    internal_stress();//only update internal stress, external stays the analytical one
    propagate();
    if (iTime%100 != 0)  continue;
    io::write_vectors("config_"+to_string(iTime)+".dat", 
                      {&bin_center, &disp, &int_stress, &ext_stress}, "r\tdisp\tint_stress\text_stress");
  }
  io::write_vectors("disp_hertz.dat", {&bin_center, &disp, &solution, &indenter::height}, "r\tdisp(sim)\tdisp(ana)");//*/
}



void init() {
  elast::init_uniform_bins(1, 1);
  elast::init_stiffness();
};


void test_interaction() {
  double dz = 1e-2*interaction::range;
  indenter::z0 = -interaction::range;

  ofstream output("testPot.dat");
  output << "# interaction::surfEnerg = " << interaction::surfEnerg << "\n";
  output << "# interaction::curvature = " << interaction::curvature << "\n";
  output << "# interaction::range = " << interaction::range << "\n";
  output << "\n# z\tpot1(z)\tpot2(z)\n";

  while (indenter::z0 < 2*interaction::range) {
    output << indenter::z0 << "\t" 
      << interaction::potential1(indenter::z0) << "\t" 
      << interaction::potential2(indenter::z0) << "\n";
    indenter::z0 += dz;
  }
  output.close();
};


void test_indenter() {
  elast::init_uniform_bins(100, 1);
  indenter::z0 = 0;
  indenter::init_height(elast::bin_center, &(indenter::polynomial));
  io::write_vectors("testIndenter.dat", {&elast::bin_center, &indenter::height}, "# r\tz(r)");
}


void test_Verlet() {
  elast::init_uniform_bins(1, 2);
  elast::init_stiffness();
  elast::stiffness_array[0][0] = 1.;
  elast::disp[0] = 1.; // u(t=0) relaxation
  elast::disp_old[0] = 0.999;

  ofstream output("verlet.dat");
  output << "# time\tdisplacement\n";
  damping = 0.0;
  dTime = 0.1;
  uint32_t nTime = 1000;
  for (int iTime = 0; iTime < nTime; ++iTime) {
    elast::set_stress(0);
    elast::internal_stress();
    elast::propagate();
    output << iTime*dTime << "\t" << elast::disp[0] << "\n";
  }
  output.close();
}





void test_Hertz() {
  elast::init_uniform_bins(250, 2.);
  elast::init_stiffness();
  elast::set_disp(0);
  elast::set_disp_old(0);
  indenter::init_height(elast::bin_center, &(indenter::polynomial));
  io::write_array("stiff.dat", elast::stiffness_array, "v sBin v || < uBin >");
  io::write_vectors("indenter.dat", {&elast::bin_center, &indenter::height}, "r\tz(r)");

  uint32_t nTime = 100000;
  for (uint32_t iTime = 0; iTime < nTime; ++iTime) {
    double time = iTime*dTime;
    //elast::set_disp(1.);
    elast::set_stress(0);
    elast::compute_stress();
    elast::propagate();
    if (iTime%1000) continue;
    io::write_vectors("disp_"+to_string(iTime)+".dat", 
      {&elast::bin_center, &elast::disp, &elast::int_stress, &elast::ext_stress}, "r\tdisp\tint_stress\text_stress");
  }
}


int main() {
  test_ellip_int();
  elast::test_stress();
  elast::test_disp_propagation();
  test_indenter();
  cout << boost::math::ellint_2(0.5) << "\n";
  return 0;
}