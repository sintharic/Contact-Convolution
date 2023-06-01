#include "Indenter.h"

extern const uint8_t UNIFORM, LOG, DYNAMIC;
extern const uint8_t FLAT, POLY, SPHERE;
extern const geometry HERTZ;

double Indenter::pos(double z) {
  return function(z, geom);
}

void Indenter::init(vector<double>& edges, vector<double>& centers, map<string, string> toml) {
  for (const auto& [key, value] : toml) {
    //cout << key << " " << value << endl;//DEBUG
    size_t idx = key.find(".");
    if (idx != string::npos && key.substr(0,idx) != "Indenter") continue;//ADJUST for ID
    string name = key.substr(idx+1);
    if (name == "type") geom.type = stoi(value);
    else if (name == "radius") geom.radius = stod(value);
    else if (name == "exponent") geom.exponent = stod(value);
    else if (name == "z0") geom.z0 = stod(value);
    //else std::cout << "# unknown geometry parameter: " << name << std::endl;//DEBUG
  }
  init(edges, centers);
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

Indenter::Indenter(vector<double>& edges, vector<double>& centers, map<string,string> toml) {
  init(edges, centers, toml);
}

void Indenter::write_params(const string& filename) {
  ofstream output(filename, std::ios::app);
  if (!output.is_open()) return;

  output << "[Indenter]\n";
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