#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

double xmax = 20;

double integrand(double x, double phi) {
  return 1./sqrt(1 - cos(phi)*2*x/(x*x + 1));
};

double intK(double x, uint64_t nphi) {
  double dphi = 2*M_PI/nphi;
  double integral = 0.;
  for (int ix = 0; ix < nphi; ++ix) {
    double phi = (ix+0.5)*dphi;
    integral += integrand(x,phi)*dphi;
  }
  return integral;
};

void computeK(uint64_t nx, uint64_t nphi, string filename) {
  ofstream file(filename);
  file << "# r/r'\tK(r/r')\tE'*G(r/r')\n";
  //uint64_t nx = 10000, nphi = 100;
  double dx = xmax/nx;
  for (int ix = 0; ix < nx; ++ix) {
    double x = (ix+0.5)*dx;
    double K = intK(x,nphi);
    file << x << "\t" << K << "\t" << K/sqrt(x*x + 1)/(2*M_PI*M_PI) << "\n";
  }
  file.close();
};

int main() {
  //computeK(1000, 100, "GF_3_2.dat");
  //computeK(10000, 100, "GF_4_2.dat");
  //computeK(100000, 100, "GF_5_2.dat");
  //computeK(10000, 1000, "GF_4_3.dat");
  //computeK(10000, 10000, "GF_4_4.dat");
  //computeK(10000, 100000, "GF_4_5.dat");
  computeK(100000, 100000, "GF_5_5.dat");

  return 0;
}