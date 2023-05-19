#include "Indenter.h"

extern const uint8_t UNIFORM, LOG, DYNAMIC;
extern const uint8_t FLAT, POLY, SPHERE;
extern const geometry HERTZ;

double Indenter::pos(double z) {
  return function(z, geom);
}

void Indenter::init(vector<double>& rval, geometry g) {
  geom = g;
  bin_center.resize(rval.size());
  height.resize(rval.size());

  if (geom.type == POLY) function = &(polynomial);
  else if (geom.type == FLAT) {
    cerr << "flat punch indenter not yet implemented.\n";
    exit(1);
  }
  else if (geom.type == SPHERE) {
    cerr << "spherical indenter not yet implemented.\n";
    exit(1);
  }

  for (int ir = 0; ir < rval.size(); ++ir) {
    bin_center[ir] = rval[ir];
    height[ir] = function(rval[ir], geom);
  }
};

Indenter::Indenter(vector<double>& rval, geometry g) {
  init(rval, g);
}