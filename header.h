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

static vector<double> uniform_bins(uint32_t N, double Rmax) {
/* 
  discretizes radial coordinate into <N> equidistant points from 0 to <Rmax>. 
*/
  vector<double> edges(N+1);

  double dr = Rmax/N;
  for (int iEdge = 0; iEdge <= N; ++iEdge) 
    edges[iEdge] = iEdge*dr;
  
  return edges;
};

static vector<double> bins(uint32_t N, double Rmax, uint8_t kind) {
  if (kind==UNIFORM) return uniform_bins(N, Rmax);
  else {
    cerr << "non-uniform grids not implemented yet.\n";
    exit(1);
  }
}

static vector<double> arithmetic_centers(vector<double>& edges) {
  vector<double> result(edges.size()-1);
  for (int iBin = 0; iBin < result.size(); ++iBin) {
    result[iBin] = 0.5*(edges[iBin+1] + edges[iBin]);
  }
  return result;
}


struct geometry {
  uint8_t type = POLY;
  double z0 = 0, radius = 1., exponent = 2.;

  geometry() {};
  geometry(map<string,string> params) {
    for (const auto& [key, value] : params) {
      cout << key << " " << value << endl;//DEBUG
      if (key == "radius") radius = stod(value);
      else if (key == "exponent") exponent = stod(value);
      else if (key == "z0") z0 = stod(value);
      else if (key == "type") type = stoi(value);
      else std::cout << "# unknown geometry parameter: " << key << std::endl;//DEBUG
    }
  }
};

struct potparams {
  double curvature = 10.;
  double surfEnerg = 1e-4;
  double range = M_PI*sqrt(surfEnerg/(2*curvature));
  double potCurveRel = -1;

  potparams() {};
  potparams(double c, double gamma) : curvature(c), surfEnerg(gamma) {
    range = M_PI*sqrt(surfEnerg/(2*curvature));
  }
  potparams(map<string,string> params) {
    for (const auto& [key, value] : params) {
      //cout << key << " " << value << endl;//DEBUG
      if (key == "curvature") curvature = stod(value);
      else if (key == "surfEnerg") surfEnerg = stod(value);
      else if (key == "range") range = stod(value);
      else if (key == "potCurveRel") potCurveRel = stod(value);
      else std::cout << "# unknown potential parameter: " << key << std::endl;//DEBUG
    }
  }
};

namespace io {
  vector< map<string,string> > read_pairs(const string&);
  template <typename T> vector<T> read_column(const string& filename, uint32_t col);
  double column_increment(const string&, uint32_t);
  void write_vectors(const string& filename, vector< vector<double>* > arrays, const string& header);
  void write_array(const string& filename, vector< vector<double> >& array, const string& header);
}

const geometry HERTZ;
const potparams DEFAULT_POTENTIAL(10., 1e-4);
static map<uint8_t, string> GEOMETRY = {{FLAT, "FLAT"}, {POLY, "POLY"}, {SPHERE, "SPHERE"}};
static map<uint8_t, string> GRID = {{UNIFORM, "UNIFORM"}, {LOG, "LOG"}, {DYNAMIC, "DYNAMIC"}};
static map<string, string> SUBS = {
  {"UNIFORM", to_string(UNIFORM)}, 
  {"LOG", to_string(LOG)}, 
  {"DYNAMIC", to_string(DYNAMIC)},
  {"FLAT", to_string(FLAT)},
  {"POLY", to_string(POLY)},
  {"SPHERE", to_string(SPHERE)}
};

#endif