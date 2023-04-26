#include "header.h"

double contMod = 1.;

vector<double> stress, rStress;
vector<double> dispNow, dispOld, rDisp;
vector<vector<double> > stiffness;

void init_stiffness() {
  vector<double> numericGF = io::read_column<double>("GF_5_5.dat", 2);
  double dx = io::column_increment("GF_5_5.dat", 0);
  
  //DEBUG
  cout << "# increment: " << dx << "\n";
  cout << "# GF data pts: " << numericGF.size() << "\n";
  //cout << "# GF head pts: " << numericGF[0] << " " << numericGF[1] << "\n";
  //cout << "# GF tail pts: " << numericGF[numericGF.size()-2] << " " << numericGF[numericGF.size()-1] << "\n";
  //for (auto val : numericGF) cout << val << "\n";

  for (uint32_t iStress=0; iStress<rStress.size(); ++iStress) {
    for (uint32_t iDisp=0; iDisp<rDisp.size(); ++iDisp) {
      double x = rDisp[iDisp]/rStress[iStress];
      if (x > 5) stiffness[iStress][iDisp] = contMod*M_PI*x;
      else if (x < 0.2) stiffness[iStress][iDisp] = contMod*M_PI/(1. + x*x/4);
      //else if (x > 0.99 && x < 1.01) stiffness[iStress][iDisp] = 0;
      else {
        uint32_t idx = x/dx; //TODO: more exact (interpolation)?
        stiffness[iStress][iDisp] = contMod/numericGF[idx];
      }
    }
  }
};

void init_GFtest() {
  uint32_t nx = 1000;
  stress.resize(nx/4);
  dispNow.resize(nx);
  rStress.resize(nx/4);
  rDisp.resize(nx);
  stiffness.resize(nx/4);
  for (uint32_t i=0; i<stiffness.size(); ++i) stiffness[i].resize(nx);

  //TEST
  double dxDisp = 0.02;
  for (int ix = 0; ix < nx/4; ++ix) rStress[ix] = ix*dxDisp;
  for (int ix = 0; ix < nx; ++ix) rDisp[ix] = (ix+0.5)*dxDisp;
}

void run_GFtest() {
  ofstream file("testGF.dat");
  for (int i = 0; i < stiffness[50].size(); ++i) {
    file << rDisp[i] << "\t" << 1./stiffness[50][i] << "\n";
  }
  file.close();
}

void stressConvolution() {
  for (int iStress = 0; iStress < rStress.size(); ++iStress) {
    //TODO: better determination of dr'
    double drStress;
    if (iStress>0 && iStress<rStress.size()-1) 
      drStress = 0.5*(rStress[iStress+1] - rStress[iStress-1]);
    else if (iStress > 0) 
      drStress = rStress[iStress+1] - rStress[iStress];
    else 
      drStress = rStress[iStress] - rStress[iStress-1];

    stress[iStress] = 0;
    for (int iDisp = 0; iDisp < rDisp.size(); ++iDisp) {
      //TODO: better determination of dr
      double drDisp;
      if (iDisp>0 && iDisp<rDisp.size()-1) 
        drDisp = 0.5*(rDisp[iDisp+1] - rDisp[iDisp-1]);
      else if (iDisp > 0) 
        drDisp = rDisp[iDisp+1] - rDisp[iDisp];
      else 
        drDisp = rDisp[iDisp] - rDisp[iDisp-1];

      //TODO: what to do with dr, dr'?
      stress[iStress] += stiffness[iStress][iDisp]*dispNow[iDisp];
    }
  }
}

int main() {
  init_GFtest();
  init_stiffness();
  run_GFtest();

  return 0;
}