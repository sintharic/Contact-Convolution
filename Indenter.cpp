#include "Indenter.h"

extern const uint8_t UNIFORM, LOG, DYNAMIC;
extern const uint8_t FLAT, POLY, SPHERE;
extern const geometry HERTZ;

double Indenter::pos(double z) {
  return function(z, geom);
}

void Indenter::init(vector<double>& edges, vector<double>& centers) {
  bin_center.resize(centers.size());
  height.resize(centers.size());

  if (geom.type == POLY) function = &(polynomial);
  else if (geom.type == FLAT) {
    cerr << "flat punch indenter not yet implemented.\n";
    exit(1);
  }
  else if (geom.type == SPHERE) {
    cerr << "spherical indenter not yet implemented.\n";
    exit(1);
  }

  for (int ir = 0; ir < centers.size(); ++ir) {
    bin_center[ir] = centers[ir];
    height[ir] = function(centers[ir], geom);
  }
};

Indenter::Indenter(vector<double>& edges, vector<double>& centers, geometry g) : geom(g) {
  init(edges, centers);
}