#ifndef INTER_H
#define INTER_H
//endif at EOF

#include "header.h"
#include "ElasticBody.h"
#include "Indenter.h"

class Interaction {
  potparams params;
  vector<double> *pos_ind;
  vector<double> *pos_elast;
  bool initialized = 0;
  double (*function)(double, potparams);

public:
  Interaction() {};
  Interaction(ElasticBody& elast, Indenter& ind, potparams);
  Interaction(ElasticBody& elast, Indenter& ind, map<string,string>);
  void init(ElasticBody& elast, Indenter& ind);
  void add_stress(vector<double>& stress);
  double potential(double gap) {return function(gap, params);};

  void write_params(const string&);

  static double potential1(double gap, potparams p) {
    if (gap > 0) return 0;
    return p.curvature*gap;
  };

  static double potential2(double gap, potparams p) {
    if (gap > p.range) return 0;
    else if (gap > 0) return p.surfEnerg*M_PI/p.range*sin(gap*M_PI/p.range)/2.;
    return p.curvature*gap;
  };

  static void test_potentials() {
    cout << "# Interaction::test_potentials()\n";//FLOW
    geometry g;
    double dz = 1e-2*DEFAULT_POTENTIAL.range;

    ofstream output("potential.dat");
    output << "# params.surfEnerg = " << DEFAULT_POTENTIAL.surfEnerg << "\n";
    output << "# params.curvature = " << DEFAULT_POTENTIAL.curvature << "\n";
    output << "# params.range = " << DEFAULT_POTENTIAL.range << "\n";
    output << "\n# z\tpot1(z)\tpot2(z)\n";

    while (g.z0 < 2*DEFAULT_POTENTIAL.range) {
      output << g.z0 << "\t" 
        << potential1(g.z0, DEFAULT_POTENTIAL) << "\t" 
        << potential2(g.z0, DEFAULT_POTENTIAL) << "\n";
      g.z0 += dz;
    }
    output.close();
  };
};

#endif