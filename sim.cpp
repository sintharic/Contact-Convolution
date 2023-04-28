#include "header.h"

double damping = 2.; // 2.0 is exactly critical damping

namespace elast {
  double contMod = 1.;
  vector<double> stress, rStress;
  vector<double> dispNow, dispOld, rDisp;
  vector<vector<double> > GreensFunction;

  void init_GF() {
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
        if (x > 5) GreensFunction[iStress][iDisp] = 1./(contMod*M_PI*x);
        else if (x < 0.2) GreensFunction[iStress][iDisp] = (1. + x*x/4)/(contMod*M_PI);
        //else if (x > 0.99 && x < 1.01) GreensFunction[iStress][iDisp] = 0;
        else {
          uint32_t idx = x/dx; //TODO: more exact (interpolation)?
          GreensFunction[iStress][iDisp] = numericGF[idx]/contMod;
        }
      }
    }
  };

  void init_stiffness () {

  };

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
        stress[iStress] += dispNow[iDisp]/GreensFunction[iStress][iDisp];
      }
    }
  }
}

namespace indenter {
  double z0 = 0;
  double dzInd = 1e-3;
  double pressure = 0.2;
  double radius = 1.;
  double exponent = 2.;
  vector<double> height;

  double get_z(double r) {
    double h = pow(r*r/radius, exponent/2)/exponent;
    return(z0 + h);
  };

  void init_height(vector<double>& rval) {
    height.resize(rval.size());
    for (int ir = 0; ir < rval.size(); ++ir) {
      height[ir] = get_z(rval[ir]);
    }
  };
}

namespace interaction {
  double curvature = 1.;
  double surfEnerg = 1e-4;
  double range = M_PI*sqrt(surfEnerg/(2*curvature));

  double potential1(double gap) {
    if (gap > 0) return 0;
    return curvature*gap;
  }

  double potential2(double gap) {
    if (gap > range) return 0;
    else if (gap > 0) return surfEnerg*M_PI/range*sin(gap*M_PI/range)/2.;
    return curvature*gap;
  }
}

void init_GFtest() {
  uint32_t nx = 1000;
  elast::stress.resize(nx/4);
  elast::dispNow.resize(nx);
  elast::rStress.resize(nx/4);
  elast::rDisp.resize(nx);
  elast::GreensFunction.resize(nx/4);
  for (uint32_t i=0; i<elast::GreensFunction.size(); ++i) elast::GreensFunction[i].resize(nx);

  //TEST
  double dxDisp = 0.02;
  for (int ix = 0; ix < nx/4; ++ix) elast::rStress[ix] = ix*dxDisp;
  for (int ix = 0; ix < nx; ++ix) elast::rDisp[ix] = (ix+0.5)*dxDisp;
}

void run_GFtest() {
  ofstream file("testGF.dat");
  for (int i = 0; i < elast::GreensFunction[50].size(); ++i) {
    file << elast::rDisp[i] << "\t" << elast::GreensFunction[50][i] << "\n";
  }
  file.close();
}

void propagate() {
  return;
}

int main() {
  init_GFtest();
  elast::init_GF();
  run_GFtest();

  return 0;
}