#include "header.h"
#include "ElasticBody.h"
#include "Indenter.h"
#include "Interaction.h"

extern const uint8_t UNIFORM, LOG, DYNAMIC;
extern const uint8_t FLAT, POLY, SPHERE;
extern const geometry HERTZ;


class Simulation {
  vector<double> bin_edge, bin_center;
  ElasticBody elastic;
  Indenter indenter;
  Interaction interaction;
  
  uint32_t nTime = 10, nBin = 100, frameInterval = 10;
  uint8_t grid = UNIFORM;
  double dTime = 0.05, Rmax = 8., area = M_PI*Rmax*Rmax;

public:
  void init(const string&);
  void write_params(const string&);
  void run(), update_bins();
};

void Simulation::init(const string& filename) {
  cerr << "# Simulation::init()\n";//FLOW

  auto maps = io::read_pairs(filename);
  map<string,string> p_global=maps[0], p_elast=maps[1], p_indenter=maps[2], p_inter=maps[3];

  for (const auto& [key, value] : p_global) {
    cout << key << " " << value << endl;//DEBUG
    if (key == "Rmax") Rmax = stod(value);
    else if (key == "dTime") dTime = stod(value);
    else if (key == "nTime") nTime = stoi(value);
    else if (key == "nBin") nBin = stoi(value);
    else if (key == "grid") grid = stoi(value);
    else if (key == "frameInterval") frameInterval = stoi(value);
    else cout << "# unknown elastic parameter: " << key << endl;//DEBUG
  }
  area = M_PI*Rmax*Rmax;

  bin_edge = bins(nBin, Rmax, grid);
  bin_center = arithmetic_centers(bin_edge);

  elastic = ElasticBody(bin_edge, bin_center, p_elast);
  indenter = Indenter(bin_edge, bin_center, p_indenter);
  interaction = Interaction(elastic, indenter, p_inter);

  write_params("params.out");
}

void Simulation::write_params(const string& filename) {
  ofstream output(filename);
  if (!output.is_open()) return;

  output << "Global:\n";
  output << "  nBin = " << nBin << "\n";
  output << "  Rmax = " << Rmax << "\n";
  output << "  grid = " << GRID[grid] << "\n";
  output << "  nTime = " << nTime << "\n";
  output << "  dTime = " << dTime << "\n";
  output << "  frameInterval = " << frameInterval << "\n";
  output << "\n";
  output.close();

  elastic.write_params(filename);
  indenter.write_params(filename);
  interaction.write_params(filename);
}

void Simulation::run() {
  indenter.write_config();
  for (int iTime = 0; iTime < nTime; ++iTime) {
    elastic.set_stress(0);
    elastic.stress_internal();
    elastic.stress_external();
    interaction.add_stress(elastic.ext_stress);
    elastic.propagate(dTime);
    if (grid == DYNAMIC) update_bins();
    if (iTime%frameInterval) continue;
    elastic.write_config(iTime);
  }
  elastic.write_config();
}

void Simulation::update_bins() {
  terminate("DYNAMIC bins not implemented yet.");
  /*uint32_t idx = interaction.idx_contact();
  bin_edge = dynamic_bins(nBin, idx, Rmax);
  bin_center = arithmetic_centers(bin_edge);
  elastic.update_bins(bin_edge, bin_center);
  indenter.update_bins(bin_edge, bin_center);//*/
}


int main() {
  Simulation sim;
  sim.init("params.in");
  sim.run();
}





// ---------------------------------





int old_main() {
  double damping = 1.8;
  double dTime = 0.05;

  ElasticBody elast(250, 8., UNIFORM);
  elast.set_damping(damping);
  elast.set_disp(0);
  elast.set_disp_old(0);
  Indenter indenter(elast.bin_edge, elast.bin_center, HERTZ);
  io::write_array("stiff_hertz.dat", elast.stiffness_array, "v sBin v || < uBin >");
  //io::write_vectors("indenter.dat", {&elast.bin_center, &indenter.height}, "r\tz(r)");
  Interaction inter(elast, indenter, potparams(elast.get_stiff(), 0));

  uint32_t nTime = 10000;
  for (uint32_t iTime = 0; iTime < nTime; ++iTime) {
    double time = iTime*dTime;
    //elast.set_disp(1.);
    elast.set_stress(0);
    elast.stress_internal();
    elast.stress_external();
    inter.add_stress(elast.ext_stress);
    elast.propagate(dTime);
    if (iTime%100) continue;
    io::write_vectors("disp_"+to_string(iTime)+".dat", 
      {&elast.bin_center, &elast.disp, &elast.int_stress, &elast.ext_stress}, "r\tdisp\tint_stress\text_stress");
  }
  return 0;
}