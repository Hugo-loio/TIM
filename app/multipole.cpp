#include <iostream>
#include <fstream>
#include <armadillo>
#include "OData.h"
#include "SSH.h"
#include "SSH2D.h"
#include "BBH2D.h"

using namespace std;
using namespace arma;

void scanQuadrupoleBBH2D(int nPoints, int nx, int ny, char* argv0, string name){
  OData out(argv0, name);
  double delta = 2/(double)nPoints;
  vector<vector<double>> data;
  vector<double> t1;
  vector<double> quad;
  BBH2D bbh2D(2,1);
  for(int i = 0; i <= nPoints; i++){
    t1.push_back(i*delta);
    bbh2D.setIntraHop(i*delta);
    quad.push_back(bbh2D.getQuadrupoleNested(nx,ny));
  }
  data.push_back(t1);
  data.push_back(quad);
  out.data(data);
}


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

  scanQuadrupoleBBH2D(100, 20, 20, argv[0], "QuadrupoleBBH2D_20_20.dat");

  return 0;
}
