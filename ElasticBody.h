#ifndef ELASTIC_H
#define ELASTIC_H
//endif at EOF

#include "header.h"

class ElasticBody {
  double Estar = 1., Estar_inf = 1.;
  double area = M_PI*4;
  uint32_t nBin = 32;
  uint32_t ID;

  // Maxwell parameters
  uint32_t nMaxwell = 0;
  vector<double> invTauMw, stiffFacMw;
  vector<vector<vector<double> > > dispMw;

  // experimental parameters (force/displacement control)
  uint8_t fConstraint = 0;
  double force = 0.1, zPos = 0.;

  // solver parameters
  double massGFMD = 1., damping = 1.8;

  void init_uniform_bins(uint32_t N, double Rmax);
  void init_fields_from_bins();
  void init_stiffness();
  void fit_stiffness();
  void read_stiffness(const string&);
  double max_stiff = 0, dispCOM = 0;

  ofstream moni;

public:
  uint8_t initialized = 0;
  vector<double> bin_edge; // length nBin+1
  vector<double> bin_center, bin_area; // length nBin
  vector<double> ext_stress, int_stress, disp, disp_old; // length nBin
  vector<vector<double> > stiffness_array;

  //new
  ElasticBody() : ID(0) {};
  ElasticBody(vector<double>&, vector<double>&, map<string,string>);
  ElasticBody(uint32_t i) : ID(i) {};
  ElasticBody(vector<double>&, vector<double>&, map<string,string>, uint32_t);
  void init(vector<double>&, vector<double>&, map<string,string>);
  //old
  ElasticBody(uint32_t, double, uint8_t);
  ElasticBody(map<string,string>, map<string,string>);
  void init(uint32_t, double, uint8_t);
  
  void init_Maxwell(double);
  void set_stress(double val);
  void set_disp(double val);
  void set_values(vector<double>, vector<double>*);
  void set_disp_old(double val);
  void set_damping(double);
  void update_bins(vector<double>&, vector<double>&);
  void stress_internal();
  void stress_external();
  void propagate(double);

  void write_config();
  void write_config(uint32_t);
  void write_params(const string&);

  // access
  double get_Estar() {return Estar;};
  double get_force() {return force;};
  double get_area() {return area;};
  double get_stiff() {return max_stiff;};

  // tests
  void test_Verlet();
  void test_stress();
  void test_disp_propagation();
  static void test_ellip_int();
};

#endif