#ifndef HEADER_H
#define HEADER_H
//endif at EOF

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>

using namespace std;

namespace io {
  template <typename T> vector<T> read_column(const string& filename, uint32_t col);
  double column_increment(const string&, uint32_t);
  void write_vectors(const string& filename, vector< vector<double>* > arrays, const string& header);
  void write_array(const string& filename, vector< vector<double> >& array, const string& header);
}

/*namespace elast {
  double Estar = 1.;
  uint32_t nBin = 32;
  vector<double> bin_center(nBin), bin_area(nBin); // length nBin
  vector<double> stress(nBin), disp(nBin), disp_old(nBin); // length nBin
  vector<double> bin_edge; // length nBin+1
  vector<vector<double> > stiffness;

  void init_uniform_pos(uint32_t N, double Rmax);
  void init_stiffness ();
  void set_stress(double val);
  void internal_stress();
  void propagate();
}//*/

#endif