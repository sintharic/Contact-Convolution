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

  for (int iBin = 0; iBin < centers.size(); ++iBin) {
    bin_center[iBin] = centers[iBin];
    height[iBin] = function(centers[iBin], geom);
  }
};

Indenter::Indenter(vector<double>& edges, vector<double>& centers, geometry g) : geom(g) {
  init(edges, centers);
}

void Indenter::write_params(const string& filename) {
  ofstream output(filename, std::ios::app);
  if (!output.is_open()) return;

  output << "Indenter:\n";
  output << "  type = " << GEOMETRY[geom.type] << "\n";
  if (geom.type == POLY) output << "  exponent = " << geom.exponent << "\n";
  output << "  radius = " << geom.radius << "\n";
  output << "  z0 = " << geom.z0 << "\n";
  output << "\n";

  output.close();
}

void Indenter::write_config() {
  io::write_vectors("indenter.dat", {&bin_center, &height}, "# r\tz(r)");
}

void Indenter::write_config(uint32_t iTime) {
  io::write_vectors("indenter."+to_string(iTime)+".dat", {&bin_center, &height}, "# r\tz(r)");
}