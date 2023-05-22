#ifndef INDENTER_H
#define INDENTER_H
//endif at EOF

#include "header.h"

class Indenter {
  geometry geom;
  double (*function)(double, geometry);

public:
  vector<double> bin_center, height;

  Indenter() {};
  Indenter(vector<double>&, vector<double>&, geometry);
  void init(vector<double>&, vector<double>&);
  void update_bins(vector<double>&, vector<double>&);
  double pos(double r);

  static double polynomial(double r, geometry geom) {
    double h = pow(r*r/geom.radius, geom.exponent/2)/geom.exponent;
    return(geom.z0 + h);
  };
};

#endif