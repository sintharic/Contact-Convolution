#include "header.h"
#include "ElasticBody.h"
#include "Indenter.h"
#include "Interaction.h"

extern const uint8_t UNIFORM, LOG, DYNAMIC;
extern const uint8_t FLAT, POLY, SPHERE;
extern const geometry HERTZ;



int main() {
  double damping = 1.8;
  double dTime = 0.05;

  ElasticBody elast(250, 8., UNIFORM);
  elast.set_damping(damping);
  elast.set_disp(0);
  elast.set_disp_old(0);
  Indenter ind(elast.bin_center, HERTZ);
  io::write_array("stiff_hertz.dat", elast.stiffness_array, "v sBin v || < uBin >");
  //io::write_vectors("indenter.dat", {&elast.bin_center, &ind.height}, "r\tz(r)");
  Interaction inter(elast, ind, potparams(elast.get_stiff(), 0));

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