#include "Interaction.h"

extern const uint8_t UNIFORM, LOG, DYNAMIC;
extern const uint8_t FLAT, POLY, SPHERE;
extern const geometry HERTZ;

void Interaction::init(ElasticBody& elast, Indenter& ind) {
  pos_ind = &(ind.height);
  pos_elast = &(elast.disp);
  
  if (params.surfEnerg == 0) function = potential1;
  else function = potential2;

  initialized = 1;
}

Interaction::Interaction(ElasticBody& elast, Indenter& ind, potparams p) : params(p) {
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