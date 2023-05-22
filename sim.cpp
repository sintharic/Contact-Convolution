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
  
  uint32_t nTime = 10, nBin = 100;
  uint8_t grid = UNIFORM;
  double dTime = 0.05, Rmax = 8., area = M_PI*Rmax*Rmax;

public:
  void init(const string&);
  void run(), update_bins();
};

void Simulation::init(const string& filename) {
  auto maps = io::read_pairs("params.test");
  map<string,string> p_global=maps[0], p_elast=maps[1], p_indenter=maps[2], p_inter=maps[3];

  if (p_global.find("Rmax") != p_global.end()) Rmax = stod(p_global["Rmax"]);
  if (p_global.find("dTime") != p_global.end()) dTime = stod(p_global["dTime"]);
  if (p_global.find("nTime") != p_global.end()) nTime = stoi(p_global["nTime"]);
  if (p_global.find("nBin") != p_global.end()) nBin = stoi(p_global["nBin"]);
  if (p_global.find("grid") != p_global.end()) grid = stoi(p_global["grid"]);
  area = M_PI*Rmax*Rmax;

  bin_edge = bins(nBin, Rmax, grid);
  bin_center = arithmetic_centers(bin_edge);

  elastic = ElasticBody(bin_edge, bin_center, p_elast);
  indenter = Indenter(bin_edge, bin_center, p_indenter);
  interaction = Interaction(elastic, indenter, p_inter);
}

void Simulation::run() {
  for (int iTime = 0; iTime < nTime; ++iTime) {
    elastic.set_stress(0);
    elastic.internal_stress();
    elastic.external_stress();
    interaction.add_stress(elastic.ext_stress);
    elastic.propagate(dTime);
    if (grid == DYNAMIC) update_bins();
    if (iTime%100) continue;
    io::write_vectors("disp_"+to_string(iTime)+".dat", 
      {&elastic.bin_center, &elastic.disp, &elastic.int_stress, &elastic.ext_stress}, "r\tdisp\tint_stress\text_stress");
  }
}

void Simulation::update_bins() {
  /*uint32_t idx = interaction.idx_contact();
  bin_edge = dynamic_bins(nBin, idx, Rmax);
  bin_center = arithmetic_centers(bin_edge);
  elastic.update_bins(bin_edge, bin_center);
  indenter.update_bins(bin_edge, bin_center);//*/
}

int main() {
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
    elast.internal_stress();
    elast.external_stress();
    inter.add_stress(elast.ext_stress);
    elast.propagate(dTime);
    if (iTime%100) continue;
    io::write_vectors("disp_"+to_string(iTime)+".dat", 
      {&elast.bin_center, &elast.disp, &elast.int_stress, &elast.ext_stress}, "r\tdisp\tint_stress\text_stress");
  }
  return 0;
}