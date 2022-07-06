#include <iostream>
#include <fstream>
#include <armadillo>
#include "OData.h"
#include "BBH2D.h"
#include "BBH3D.h"
#include "DisorderedBBH3D.h"
#include "MultiThread.h"
#include <chrono>

using namespace std;
using namespace arma;

void threadPolBBH2D(int dir, int * l, vector<double> & res, vector<double> params){
  BBH2D mod(0,1);
  double pol;
  for(int i = 0; i < params.size(); i++){
    try{
      mod.setIntraHop(params[i]);
      pol = mod.getBoundPolarization(l,dir);
      res.push_back(params[i]);
      res.push_back(pol);
    }
    catch(const runtime_error & error){
      cout << "Singular matrix found at intracell hop = "  << params[i] << endl;
    }
  }
}

void threadPolBBH3D(int dir, int * l, vector<double> & res, vector<double> params){
  BBH3D mod(0,1);
  double pol;
  for(int i = 0; i < params.size(); i++){
    try{
      mod.setIntraHop(params[i]);
      pol = mod.getBoundPolarization(l,dir);
      res.push_back(params[i]);
      res.push_back(pol);
    }
    catch(const runtime_error & error){
      cout << "Singular matrix found at intracell hop = "  << params[i] << endl;
    }
  }
}

void threadQuadBBH3D(int dir, int * l, vector<double> & res, vector<double> params){
  BBH3D mod(0,1);
  double pol;
  for(int i = 0; i < params.size(); i++){
    try{
      mod.setIntraHop(params[i]);
      pol = mod.getBoundQuadrupole(l,dir);
      res.push_back(params[i]);
      res.push_back(pol);
    }
    catch(const runtime_error & error){
      cout << "Singular matrix found at intracell hop = "  << params[i] << endl;
    }
  }
}

void bbh2x1(vector<double> & res, vector<double> params){
  int l[2] = {50,50};
  threadPolBBH2D(0, l, res, params);
}

void bbh2y1(vector<double> & res, vector<double> params){
  int l[2] = {50,50};
  threadPolBBH2D(1, l, res, params);
}

void bbh3polx1(vector<double> & res, vector<double> params){
  int l[3] = {10,10,10};
  threadPolBBH3D(0, l, res, params);
}

void bbh3poly1(vector<double> & res, vector<double> params){
  int l[3] = {10,10,10};
  threadPolBBH3D(1, l, res, params);
}

void bbh3polz1(vector<double> & res, vector<double> params){
  int l[3] = {10,10,10};
  threadPolBBH3D(2, l, res, params);
}

void bbh3quadx1(vector<double> & res, vector<double> params){
  int l[3] = {10,10,10};
  threadQuadBBH3D(0, l, res, params);
}

void bbh3quady1(vector<double> & res, vector<double> params){
  int l[3] = {10,10,10};
  threadQuadBBH3D(1, l, res, params);
}

void bbh3quadz1(vector<double> & res, vector<double> params){
  int l[3] = {10,10,10};
  threadQuadBBH3D(2, l, res, params);
}

void threadQuadDisBBH3D(int dir, int * l, vector<double> & res, vector<double> params){
  DisorderedBBH3D bbh3d;
  for(int i = 0; i < params.size(); i++){
    bbh3d.setIntraHop(params[i]);
    bbh3d.setSize(l);
    bbh3d.setW(0);
    bbh3d.generateDisorder();
    res.push_back(params[i]);
    res.push_back(bbh3d.getBoundQuadrupole(dir));
  }
}

void quadDisBBH3D1(vector<double> & res, vector<double> params){
  int l[3] = {7,7,7};
  threadQuadDisBBH3D(0, l, res, params);
}

void quadDisBBH3D2(vector<double> & res, vector<double> params){
  int l[3] = {10,10,10};
  threadQuadDisBBH3D(0, l, res, params);
}


int main (int arc, char ** argv) {
  cout << "\n2D BBH" << endl;
  /*
  int l[2] = {50,100};
  BBH2D bbh2D(0.5,1);
  auto t1 = chrono::high_resolution_clock::now();
  cout << bbh2D.getBoundPolarization(l,0) << endl;
  auto t2 = chrono::high_resolution_clock::now();
  auto d1 = chrono::duration_cast<chrono::microseconds>(t2 - t1);
  cout << "Total: " << d1.count() << endl; 
  */
  //scanPolBBH2D(100, 0, l, argv[0], "BoundaryPxBBH2D_50x100.dat");
  //scanPolBBH2D(100, 1, l, argv[0], "BoundaryPyBBH2D_50x100.dat");

  vector<vector<double>> paramList;
  int nPoints = 100;
  for(int i = 0; i < nPoints; i++){
    vector<double> param; 
    param.push_back((double)i*2/(double)nPoints);
    paramList.push_back(param);
  }
  /*
  MultiThread r2pol1(bbh2x1, paramList, 8);
  r2pol1.setFile(argv[0], "BoundaryPxBBH2D_50x50.dat");
  r2pol1.run();

  MultiThread r2pol2(bbh2y1, paramList, 8);
  r2pol2.setFile(argv[0], "BoundaryPyBBH2D_50x50.dat");
  r2pol2.run();
  */

  cout << "\n3D BBH" << endl;

  /*
  MultiThread r3pol1(bbh3polx1, paramList, 8);
  r3pol1.setFile(argv[0], "BoundaryPxBBH3D_10x10x10.dat");
  r3pol1.run();

  MultiThread r3pol2(bbh3poly1, paramList, 8);
  r3pol2.setFile(argv[0], "BoundaryPyBBH3D_10x10x10.dat");
  r3pol2.run();

  MultiThread r3pol3(bbh3polz1, paramList, 8);
  r3pol3.setFile(argv[0], "BoundaryPzBBH3D_10x10x10.dat");
  r3pol3.run();

  MultiThread r3quad1(bbh3quadx1, paramList, 8);
  r3quad1.setFile(argv[0], "BoundaryQyzBBH3D_10x10x10.dat");
  r3quad1.run();

  MultiThread r3quad2(bbh3quady1, paramList, 8);
  r3quad2.setFile(argv[0], "BoundaryQxzBBH3D_10x10x10.dat");
  r3quad2.run();

  MultiThread r3quad3(bbh3quadz1, paramList, 8);
  r3quad3.setFile(argv[0], "BoundaryQxyBBH3D_10x10x10.dat");
  r3quad3.run();
  */

  //Disordered BBH3D
  cout << "\nDisorederedBBH3D" << endl;

  MultiThread rDB3D1(quadDisBBH3D1, paramList, 8);
  rDB3D1.setFile(argv[0], "BoundaryQyzNoDisBBH3D_L7.dat");
  rDB3D1.run();

  MultiThread rDB3D2(quadDisBBH3D2, paramList, 8);
  rDB3D2.setFile(argv[0], "BoundaryQyzNoDisBBH3D_L10.dat");
  rDB3D2.run();

  return 0;
}
