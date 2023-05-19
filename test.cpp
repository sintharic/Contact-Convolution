#include "header.h"
#include "ElasticBody.h"
#include "Indenter.h"
#include "Interaction.h"

extern const uint8_t UNIFORM, LOG, DYNAMIC;
extern const uint8_t FLAT, POLY, SPHERE;
extern const geometry HERTZ;

// All DOFs are set to a circular eigenfrequency of sqrt(k/m) = 1.
// Therefore, 2.0 is exactly critical damping
//double damping = 1.8;
//double dTime = 0.05;


void test_indenter() {
  cerr << "# test_indenter()\n";//FLOW
  ElasticBody elast(100, 1, UNIFORM);
  Indenter hertz(elast.bin_center, HERTZ);
  io::write_vectors("indenter.dat", {&hertz.bin_center, &hertz.height}, "# r\tz(r)");
}


void test_Verlet() {
  cerr << "# test_Verlet()\n";//FLOW
  ElasticBody elast(1, 2, UNIFORM);
  elast.stiffness_array[0][0] = 1.;
  elast.disp[0] = 1.; // u(t=0) relaxation
  elast.disp_old[0] = 0.999;

  ofstream output("verlet.dat");
  output << "# time\tdisplacement\n";
  elast.set_damping(0.0);
  double dTime = 0.1;
  uint32_t nTime = 1000;
  for (int iTime = 0; iTime < nTime; ++iTime) {
    elast.set_stress(0);
    elast.internal_stress();
    elast.propagate(dTime);
    output << iTime*dTime << "\t" << elast.disp[0] << "\n";
  }
  output.close();
}

void test_stress() {
  cerr << "# test_stress()\n";//FLOW
  ElasticBody elast(200, 8., UNIFORM);
  //elast.read_stiffness("stiff.dat");
  uint32_t nBin = elast.bin_center.size();
  vector<double> solution(nBin);

  // flat punch: calculate stress from displacement profile
  double a = 1., dist = elast.get_force()/(2.*elast.get_Estar()*a);
  for (int iBin = 0; iBin < nBin ; ++iBin) {
    if (elast.bin_center[iBin] < a) {
      elast.disp[iBin] = 0;
      solution[iBin] = (elast.get_Estar()*dist/M_PI)*pow(a*a - pow(elast.bin_center[iBin],2), -0.5);
    }
    else {
      elast.disp[iBin] = dist - (2*dist/M_PI)*asin(a/elast.bin_center[iBin]);
      solution[iBin] = 0.;
    }
  }
  elast.internal_stress();
  for (double& s : elast.int_stress) s += elast.get_force()/(elast.get_area());

  io::write_vectors("stress_flat.dat", {&elast.bin_center, &elast.int_stress, &solution, &elast.disp}, 
                    "r\tstress(sim)\tstress(ana)\tdisp(sim)");
  
  // Hertz: calculate stress from displacement profile
  double R = 1., ac = pow(0.75*R*elast.get_force()/elast.get_Estar(), 1./3);
  dist = ac*ac/R;
  for (int iBin = 0; iBin < nBin ; ++iBin) {
    if (elast.bin_center[iBin] < ac) {
      elast.disp[iBin] = pow(elast.bin_center[iBin], 2)/(2*R);
      solution[iBin] = (2.*elast.get_Estar()/M_PI/R)*sqrt(ac*ac - pow(elast.bin_center[iBin], 2));
    }
    else {
      elast.disp[iBin] = dist - (ac*ac/M_PI/R)*((2-pow(elast.bin_center[iBin]/ac,2))*asin(ac/elast.bin_center[iBin]) + sqrt(pow(elast.bin_center[iBin],2) - ac*ac)/ac);
      solution[iBin] = 0.;
    }
  }
  elast.internal_stress();
  for (double& s : elast.int_stress) s += elast.get_force()/(elast.get_area());

  io::write_vectors("config.dat", 
                        {&elast.bin_center, &elast.disp, &elast.int_stress}, "r\tdisp\tint_stress");
  io::write_vectors("stress_hertz.dat", {&elast.bin_center, &elast.int_stress, &solution, &elast.disp}, 
                    "r\tstress(sim)\tstress(ana)\tdisp(sim)");
};


void test_disp_propagation() {
  cerr << "# test_disp_propagation()\n";//FLOW
  ElasticBody elast(200, 8., UNIFORM);
  io::write_array("stiff_disp.dat", elast.stiffness_array, "v sBin v || < uBin >");
  //read_stiffness("stiff.dat");
  //io::write_array("stiff_read.dat", stiffness_array, "v sBin v || < uBin >");//DEBUG
  Indenter hertz(elast.bin_center, HERTZ);
  uint32_t nBin = elast.bin_center.size();
  vector<double> solution(nBin);

  // Hertz: calculate displacement from stress profile
  double R = 1., ac = pow(0.75*R*elast.get_force()/elast.get_Estar(), 1./3);
  double dist = ac*ac/R, pressure = elast.get_force()/elast.get_area();
  for (int iBin = 0; iBin < nBin ; ++iBin) {
    if (elast.bin_center[iBin] < ac) {
      elast.ext_stress[iBin] = 1.*pressure - (2./M_PI/R)*sqrt(ac*ac - pow(elast.bin_center[iBin], 2));
      solution[iBin] = pow(elast.bin_center[iBin], 2)/(2*R);
    }
    else {
      elast.ext_stress[iBin] = 1.*pressure;
      solution[iBin] = dist - (ac*ac/M_PI/R)*((2-pow(elast.bin_center[iBin], 2)/(ac*ac))*asin(ac/elast.bin_center[iBin]) + sqrt(pow(elast.bin_center[iBin],2) -ac*ac)/ac);
    }
  }
  uint32_t nTime = 1000;
  double dTime = 0.05;
  for (int iTime = 0; iTime < nTime; ++iTime) {
    elast.internal_stress();//only update internal stress, external stays the analytical one
    elast.propagate(dTime);
    if (iTime%10 != 0)  continue;
    io::write_vectors("config_"+to_string(iTime)+".dat", 
                      {&elast.bin_center, &elast.disp, &elast.int_stress, &elast.ext_stress}, "r\tdisp\tint_stress\text_stress");
  }
  io::write_vectors("disp_hertz.dat", {&elast.bin_center, &elast.disp, &solution, &hertz.height}, "r\tdisp(sim)\tdisp(ana)");//*/
}

void test_reader() {
  cerr << "# test_reader()\n";//FLOW
  auto maps = io::read_pairs("params.test");
  ofstream output("params_read.test");
  for (auto m : maps) {
    output << "Map:\n";
    for (const auto& [key, value] : m)
      output << key << " = " << value << "\n";
    output << "\n";
  }
  output.close();

  map<string,string> global=maps[0], elast=maps[1], indenter=maps[2], inter=maps[3];

  ElasticBody body(maps[0], maps[1]);
}

int main() {
  ElasticBody::test_ellip_int();
  Interaction::test_potentials();
  test_indenter();
  test_Verlet();
  test_stress();
  test_disp_propagation();
  test_reader();
  return 0;
}