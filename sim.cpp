#include "header.h"

// All DOFs are set to a circular eigenfrequency of sqrt(k/m) = 1.
// Therefore, 2.0 is exactly critical damping
double damping = 2.0;
double dTime = 0.01;


namespace elast {
  double Estar = 1.;
  uint32_t nBin = 32;
  vector<double> bin_center(nBin), bin_area(nBin); // length nBin
  vector<double> stress(nBin), disp(nBin), disp_old(nBin); // length nBin
  vector<double> bin_edge; // length nBin+1
  vector<vector<double> > stiffness;

  void init_uniform_bins(uint32_t N, double Rmax);
  void init_fields();
  void init_stiffness();
  void set_stress(double val);
  void internal_stress();
  void propagate();
}//*/

void elast::init_uniform_bins(uint32_t N, double Rmax) {
  /* 
  discretizes radial coordinate into <N> equidistant points from 0 to <Rmax>. 
  */

  nBin = N;
  bin_edge.resize(nBin+1);

  double dr = Rmax/nBin;
  for (int iBin = 0; iBin <= nBin; ++iBin) 
    bin_edge[iBin] = iBin*dr;
  
  init_fields();
};

void elast::init_fields() {
  /* 
  sets up all other fields according to the generated bin edges
  */

  bin_center.resize(nBin);
  bin_area.resize(nBin);
  stress.resize(nBin);
  disp.resize(nBin);
  disp_old.resize(nBin);
  stiffness.resize(nBin);

  for (int iBin = 0; iBin < nBin; ++iBin) {
    stiffness[iBin].resize(nBin);
    bin_center[iBin] = 0.5*(bin_edge[iBin+1] + bin_edge[iBin]);
    bin_area[iBin] = 2*M_PI*bin_center[iBin]*(bin_edge[iBin+1] - bin_edge[iBin]);
  }
};

void elast::init_stiffness () {
  /* 
  stiffnes represents distribution of stress per surface displacement, where 
  only a single circle of width <bin_width> at r=<bin_center> is displaced.
  */

  for (int uBin = 0; uBin < nBin; ++uBin) {
    
    // total force acting outside the displaced circle
    double Fout_per_u = 0;

    double bin_width = bin_area[uBin]/(2*M_PI*bin_center[uBin]);
    double prefac = Estar*bin_width/(2*M_PI*bin_center[uBin]*bin_center[uBin]);

    for (int sBin = 0; sBin < nBin; ++sBin) {
      if (uBin==sBin) continue;

      // equations fitted to GFMD data
      double x = bin_center[sBin]/bin_center[sBin];
      if (x>1) stiffness[sBin][uBin] = prefac*(1./((x-1)*(x-1)));
      else stiffness[sBin][uBin] = prefac*(1./((x-1)*(x-1)) + 1./((x+1)*(x+1)) + M_PI/2);
      Fout_per_u += stiffness[sBin][uBin]*bin_area[uBin];
    }

    // force inside the circle must balance out the force outside of it
    stiffness[uBin][uBin] = -Fout_per_u/bin_area[uBin];
  }
};

void elast::set_stress(double val) {
  for (int i = 0; i < stress.size(); ++i) stress[i] = val;
};

void elast::internal_stress() {
  for (int sBin = 0; sBin < stress.size(); ++sBin) {
    for (int uBin = 0; uBin < disp.size(); ++uBin) {
      stress[sBin] += disp[uBin]*stiffness[sBin][uBin];
    }
  }
};

void elast::propagate() {
  return;
};



namespace indenter {
  double z0 = 0;
  double dzInd = 1e-3;
  double pressure = 0.2;
  double radius = 1.;
  double exponent = 2.;
  vector<double> height;

  double Z(double r) {
    double h = pow(r*r/radius, exponent/2)/exponent;
    return(z0 + h);
  };

  void init_height(vector<double>& rval) {
    height.resize(rval.size());
    for (int ir = 0; ir < rval.size(); ++ir) {
      height[ir] = Z(rval[ir]);
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
  };

  double potential2(double gap) {
    if (gap > range) return 0;
    else if (gap > 0) return surfEnerg*M_PI/range*sin(gap*M_PI/range)/2.;
    return curvature*gap;
  };
}



void init() {
  elast::init_uniform_bins(1, 1);
  elast::init_stiffness();
};



void test_interaction() {
  double dz = 1e-2*interaction::range;
  indenter::z0 = -interaction::range;

  ofstream output("testPot.dat");
  output << "# interaction::surfEnerg = " << interaction::surfEnerg << "\n";
  output << "# interaction::curvature = " << interaction::curvature << "\n";
  output << "# interaction::range = " << interaction::range << "\n";
  output << "\n# z\tpot1(z)\tpot2(z)\n";

  while (indenter::z0 < 2*interaction::range) {
    output << indenter::z0 << "\t" 
      << interaction::potential1(indenter::z0) << "\t" 
      << interaction::potential2(indenter::z0) << "\n";
    indenter::z0 += dz;
  }
  output.close();
};

void test_indenter() {
  elast::init_uniform_bins(100, 1);
  indenter::z0 = 0;
  indenter::init_height(elast::bin_center);
  io::write_vectors("testIndenter.dat", {&elast::bin_center, &indenter::height}, "# r\tz(r)");
}



int main() {
  test_interaction();
  test_indenter();
  return 0;
}