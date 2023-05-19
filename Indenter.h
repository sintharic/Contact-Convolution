#ifndef INDENTER_H
#define INDENTER_H
//endif at EOF

#include "header.h"

class Indenter {
  geometry geom;
  double (*function)(double, geometry);

public:
  vector<double> bin_center, height;

  Indenter(vector<double>& rval, geometry g);
  void init(vector<double>& rval, geometry g);
  void update_bins(vector<double>& rval);
  double pos(double r);

  static double polynomial(double r, geometry geom) {
    double h = pow(r*r/geom.radius, geom.exponent/2)/geom.exponent;
    return(geom.z0 + h);
  };
};

#endif