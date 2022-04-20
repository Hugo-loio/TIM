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
  //scanQuadrupoleBBH2D(100, 20, 20, argv[0], "QuadrupoleBBH2D_20_20.dat");
  //scanQuadrupoleBBH2D(100, 10, 10, argv[0], "QuadrupoleBBH2D_10_10.dat");

  cout << "Supercell" << endl;
  BBH2D bbh2D(0.5,1);
  int l[2] = {5,5};
  int n[2] = {10,10};
  cout << "Gamma = 0.5" << endl;
  cout << "Quadrupole: " << bbh2D.getQuadrupoleNestedSupercell(l,n) << endl;
  bbh2D.setIntraHop(2);
  cout << "Gamma = 2" << endl;
  cout << "Quadrupole: " << bbh2D.getQuadrupoleNestedSupercell(l,n) << endl;

  return 0;
}
