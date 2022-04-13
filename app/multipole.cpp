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
  double k1[2] = {0, M_PI/2};
  double k2[2] = {M_PI/2, 0};
  BBH2D bbh2D(1,2);
  cout << "Intercell hop: " << 2 << endl;
  cout << "Quadrupole moment nested: " << bbh2D.getQuadrupoleNested(20,20) << endl;
  cout << "Quadrupole moment nested: " << bbh2D.getQuadrupoleNested(20,10) << endl;
  cout << "Quadrupole moment nested: " << bbh2D.getQuadrupoleNested(10,10) << endl;
  bbh2D.setInterHop(3);
  cout << "Intercell hop: " << 3 << endl;
  cout << "Quadrupole moment nested: " << bbh2D.getQuadrupoleNested(30,20) << endl;
  bbh2D.setInterHop(1.5);
  cout << "Intercell hop: " << 1.5 << endl;
  cout << "Quadrupole moment nested: " << bbh2D.getQuadrupoleNested(20,20) << endl;
  bbh2D.setInterHop(1);
  cout << "Intercell hop: " << 1 << endl;
  cout << "Quadrupole moment nested: " << bbh2D.getQuadrupoleNested(10,10) << endl;
  bbh2D.setInterHop(0.5);
  cout << "Intercell hop: " << 0.5 << endl;
  cout << "Quadrupole moment nested: " << bbh2D.getQuadrupoleNested(20,20) << endl;
  cout << "Quadrupole moment nested: " << bbh2D.getQuadrupoleNested(10,20) << endl;
  cout << "Quadrupole moment nested: " << bbh2D.getQuadrupoleNested(20,10) << endl;
  cout << "Quadrupole moment nested: " << bbh2D.getQuadrupoleNested(20,20,k1) << endl;
  cout << "Quadrupole moment nested: " << bbh2D.getQuadrupoleNested(20,20,k2) << endl;
  bbh2D.setInterHop(0.2);
  cout << "Intercell hop: " << 0.2 << endl;
  cout << "Quadrupole moment nested: " << bbh2D.getQuadrupoleNested(30,30) << endl;

  return 0;
}
