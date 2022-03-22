#include <iostream>
#include <fstream>
#include <armadillo>
#include "OData.h"
#include "SSH.h"
#include "SSH2D.h"
#include "BBH2D.h"

using namespace std;
using namespace arma;

int main (int arc, char ** argv) {
  cout << "\n2D BBH" << endl;
  BBH2D bbh2D(1,2);
  cout << "Intercell hop: " << 2 << endl;
  cout << "Quadrupole moment nested: " << bbh2D.getQuadrupoleNested(10,10) << endl;
  bbh2D.setInterHop(3);
  cout << "Intercell hop: " << 3 << endl;
  cout << "Quadrupole moment nested: " << bbh2D.getQuadrupoleNested(10,10) << endl;
  bbh2D.setInterHop(1.5);
  cout << "Intercell hop: " << 1.5 << endl;
  cout << "Quadrupole moment nested: " << bbh2D.getQuadrupoleNested(10,10) << endl;
  bbh2D.setInterHop(1);
  cout << "Intercell hop: " << 1 << endl;
  cout << "Quadrupole moment nested: " << bbh2D.getQuadrupoleNested(10,10) << endl;
  bbh2D.setInterHop(0.5);
  double k[2] = {0, -M_PI/2};
  cout << "Intercell hop: " << 0.5 << endl;
  cout << "Quadrupole moment nested: " << bbh2D.getQuadrupoleNested(20,20,k) << endl;
  bbh2D.setInterHop(0.2);
  cout << "Intercell hop: " << 0.2 << endl;
  cout << "Quadrupole moment nested: " << bbh2D.getQuadrupoleNested(10,10) << endl;

  return 0;
}
