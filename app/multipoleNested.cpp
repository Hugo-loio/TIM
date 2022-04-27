#include <iostream>
#include <fstream>
#include <armadillo>
#include "OData.h"
#include "SSH.h"
#include "SSH2D.h"
#include "BBH2D.h"
#include "BBH3D.h"

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

void scanQuadrupoleBBH2DSupercell(int nPoints, int * l, int * n, char* argv0, string name){
  OData out(argv0, name);
  double delta = 2/(double)nPoints;
  vector<vector<double>> data;
  vector<double> t1;
  vector<double> quad;
  BBH2D bbh2D(2,1);
  for(int i = 0; i <= nPoints; i++){
    t1.push_back(i*delta);
    bbh2D.setIntraHop(i*delta);
    quad.push_back(bbh2D.getQuadrupoleNestedSupercell(l,n));
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
  int l2[2] = {4,4};
  int n2[2] = {10,10};
  //scanQuadrupoleBBH2DSupercell(100, l2, n2, argv[0], "QuadrupoleBBH2DSupercell_4_4.dat");

  cout << "3D BBH" << endl;
  BBH3D bbh3D(2,1);
  int l3[3] = {2,2,2};
  int n3[3] = {5,5,5};
  cout << "Intracell 2: " << endl;
  cout << "Octupole k: " << bbh3D.getOctupoleNested(4, 4, 4) << endl;
  cout << "Octupole k: " << bbh3D.getOctupoleNested(5, 5, 5) << endl;
  cout << "Octupole supercell: " << bbh3D.getOctupoleNestedSupercell(l3, n3) << endl;
  bbh3D.setIntraHop(0.5);
  cout << "Intracell 0.5: " << endl;
  cout << "Octupole k: " << bbh3D.getOctupoleNested(4, 4, 4) << endl;
  cout << "Octupole supercell: " << bbh3D.getOctupoleNestedSupercell(l3, n3) << endl;

  return 0;
}
