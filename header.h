#ifndef HEADER_H
#define HEADER_H
//endif at EOF

#include <iostream>
#include <fstream>
#include <cmath>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/jacobi_elliptic.hpp>
#include <gsl/gsl_sf_ellint.h>
#include <string>
#include <vector>

using namespace std;

/*// All DOFs are set to a circular eigenfrequency of sqrt(k/m) = 1.
// Therefore, 2.0 is exactly critical damping
double damping = 1.8;
double dTime = 0.05;//*/

namespace io {
  template <typename T> vector<T> read_column(const string& filename, uint32_t col);
  double column_increment(const string&, uint32_t);
  void write_vectors(const string& filename, vector< vector<double>* > arrays, const string& header);
  void write_array(const string& filename, vector< vector<double> >& array, const string& header);
}

/*namespace elast {
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

/*namespace indenter {
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

/*namespace interaction {
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

#endif