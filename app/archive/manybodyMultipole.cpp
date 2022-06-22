#include <iostream>
#include <iomanip>
#include "SSH.h"
#include "BBH2D.h"
#include "BBH3D.h"
#include "SOTAI.h"
#include "OData.h"

void scanQuadrupoleBBH2D(int nPoints, int * l, char* argv0, string name){
  OData out(argv0, name);
  double delta = 2/(double)nPoints;
  vector<vector<double>> data;
  vector<double> t1;
  vector<double> quad;
  BBH2D bbh2D(2,1);
  for(int i = 0; i <= nPoints; i++){
    t1.push_back(i*delta);
    bbh2D.setIntraHop(i*delta);
    quad.push_back(bbh2D.getQuadrupoleManyBody(l));
  }
  data.push_back(t1);
  data.push_back(quad);
  out.data(data);
}

void scanOctupoleBBH3D(int nPoints, int * l, char* argv0, string name){
  OData out(argv0, name);
  double delta = 2/(double)nPoints;
  vector<vector<double>> data;
  vector<double> t1;
  vector<double> oct;
  BBH3D bbh3D(2,1);
  for(int i = 0; i <= nPoints; i++){
    t1.push_back(i*delta);
    bbh3D.setIntraHop(i*delta);
    oct.push_back(bbh3D.getOctupoleManyBody(l));
  }
  data.push_back(t1);
  data.push_back(oct);
  out.data(data);
}

void scanQuadrupoleSOTAI(int nPoints, int * l, char* argv0, string name){
  OData out(argv0, name);
  double delta = 2/(double)nPoints;
  vector<vector<double>> data;
  vector<double> m;
  vector<double> quad;
  SOTAI sotai(2);
  for(int i = 0; i <= nPoints; i++){
    m.push_back(i*delta);
    sotai.setM(i*delta);
    quad.push_back(sotai.getQuadrupoleManyBody(l));
  }
  data.push_back(m);
  data.push_back(quad);
  out.data(data);
}


int main (int arc, char ** argv) {
  //SSH
  SSH ssh(2,1);
  cout << "SSH" << endl;
  cout << "Trivial" << endl;
  cout << "Many-body: " << ssh.polarization(10) << " Berry phase: " <<  ssh.berryPhase(10) << endl;

  cout << "Topological" << endl;
  ssh.setIntraHop(0.5);
  cout << "Many-body: " << ssh.polarization(10) << " Berry phase: " <<  ssh.berryPhase(10) << endl;

  //BBH2D
  BBH2D bbh2D(2,1);
  int l2D[2] = {10,10};
  int l2D2[2] = {5,5};
  int l2D3[2] = {15,15};
  int l2D4[2] = {20,20};
  cout << "\nBBH2D" << endl;
  cout << "Trivial" << endl;
  cout << "Many-body: " << bbh2D.getQuadrupoleManyBody(l2D) << " Nested Wilson: " <<  bbh2D.getQuadrupoleNested(10,10) << endl;

  cout << "Topological" << endl;
  bbh2D.setIntraHop(0.5);
  cout << "Many-body: " << bbh2D.getQuadrupoleManyBody(l2D) << " Nested Wilson: " <<  bbh2D.getQuadrupoleNested(10,10) << endl;
  scanQuadrupoleBBH2D(100, l2D, argv[0], "QuadrupoleManyBodyBBH2D_10x10.dat");
  //scanQuadrupoleBBH2D(100, l2D2, argv[0], "QuadrupoleManyBodyBBH2D_5x5.dat");
  //scanQuadrupoleBBH2D(100, l2D3, argv[0], "QuadrupoleManyBodyBBH2D_15x15.dat");
  //scanQuadrupoleBBH2D(100, l2D4, argv[0], "QuadrupoleManyBodyBBH2D_20x20.dat");

  //BBH3D
  BBH3D bbh3D(2,1);
  int l3D[3] = {4,4,4};
  int l3D2[3] = {5,5,5};
  int l3D3[3] = {6,6,6};
  cout << "\nBBH3D" << endl;
  cout << "Trivial" << endl;
  cout << "Many-body: " << bbh3D.getOctupoleManyBody(l3D) << " Nested Wilson: " <<  bbh3D.getOctupoleNested(4,4,4) << endl;

  cout << "Topological" << endl;
  bbh3D.setIntraHop(0.5);
  cout << "Many-body: " << bbh3D.getOctupoleManyBody(l3D) << " Nested Wilson: " <<  bbh3D.getOctupoleNested(4,4,4) << endl;
  scanOctupoleBBH3D(100, l3D, argv[0], "OctupoleManyBodyBBH3D_4x4x4.dat");
  //scanOctupoleBBH3D(100, l3D2, argv[0], "OctupoleManyBodyBBH3D_5x5x5.dat");
  //scanOctupoleBBH3D(100, l3D3, argv[0], "OctupoleManyBodyBBH3D_6x6x6.dat");

  //SOTAI
  cout << "\nSOTAI" << endl;
  scanQuadrupoleSOTAI(100, l2D, argv[0], "QuadrupoleManyBodySOTAI_10x10.dat");
  //scanQuadrupoleSOTAI(100, l2D4, argv[0], "QuadrupoleManyBodySOTAI_20x20.dat");
}
