#include "Interaction.h"

extern const uint8_t UNIFORM, LOG, DYNAMIC;
extern const uint8_t FLAT, POLY, SPHERE;
extern const geometry HERTZ;

void Interaction::init(ElasticBody& elast, Indenter& ind) {
  pos_ind = &(ind.height);
  pos_elast = &(elast.disp);

  if (params.potCurveRel > 0) params.curvature = params.potCurveRel*elast.get_stiff();
  
  if (params.surfEnerg == 0) function = potential1;
  else function = potential2;

  initialized = 1;
}

Interaction::Interaction(ElasticBody& elast, Indenter& ind, potparams p) : params(p) {
  init(elast, ind);
}

Interaction::Interaction(ElasticBody& elast, Indenter& ind, map<string,string> toml) {
  for (const auto& [key, value] : toml) {
    //cout << key << " " << value << endl;//DEBUG
    size_t idx = key.find(".");
    if (idx != string::npos && key.substr(0,idx) != "Interaction") continue;//ADJUST for ID
    string name = key.substr(idx+1);
    if (name == "curvature") params.curvature = stod(value);
    else if (name == "surfEnerg") params.surfEnerg = stod(value);
    else if (name == "range") params.range = stod(value);
    else if (name == "potCurveRel") params.potCurveRel = stod(value);
    //else std::cout << "# unknown potential parameter: " << name << std::endl;//DEBUG
  }
  init(elast, ind);
}


void Interaction::add_stress(vector<double>& stress) {
  //cerr << "add_stress()\n";//FLOW
  //double pressure = force/area;
  for (int iBin = 0; iBin < stress.size(); ++iBin) {
    stress[iBin] = stress[iBin] + potential((*pos_ind)[iBin]-(*pos_elast)[iBin]);
    //ext_stress[iBin] += 1.*pressure;
  }
};


void Interaction::write_params(const string& filename) {
  ofstream output(filename, std::ios::app);
  if (!output.is_open()) return;

  output << "[Interaction]\n";
  output << "  surfEnerg = " << params.surfEnerg << "\n";
  if (params.potCurveRel > 0) {
    output << "  potCurveRel = " << params.potCurveRel << "\n";
    output << "# curvature = " << params.curvature << "\n";
    output << "# range = " << params.range << "\n";
  }
  else {
    output << "  curvature = " << params.curvature << "\n";
    output << "  range = " << params.range << "\n";
  }
  output << "\n";

  output.close();
}