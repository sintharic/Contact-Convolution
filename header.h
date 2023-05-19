#ifndef HEADER_H
#define HEADER_H
//endif at EOF

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <cmath>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/jacobi_elliptic.hpp>
#include <algorithm>

using namespace std;

static void terminate(const string& msg) {
  cerr << msg << "\n";
  exit(1);
}

const uint8_t UNIFORM = 0, LOG = 1, DYNAMIC = 2;
const uint8_t FLAT = 0, POLY = 1, SPHERE = 2;


struct geometry {
  uint8_t type = POLY;
  double z0 = 0, radius = 1., exponent = 2.;
};

struct potparams {
  double curvature = 10.;
  double surfEnerg = 1e-4;
  double range = M_PI*sqrt(surfEnerg/(2*curvature));
  potparams();
  potparams(double c, double gamma) : curvature(c), surfEnerg(gamma) {
    range = M_PI*sqrt(surfEnerg/(2*curvature));
  }
};

const geometry HERTZ;
const potparams DEFAULT_POTENTIAL(10., 1e-4);

namespace io {
  vector< map<string,string> > read_pairs(const string&);
  template <typename T> vector<T> read_column(const string& filename, uint32_t col);
  double column_increment(const string&, uint32_t);
  void write_vectors(const string& filename, vector< vector<double>* > arrays, const string& header);
  void write_array(const string& filename, vector< vector<double> >& array, const string& header);
}

#endif