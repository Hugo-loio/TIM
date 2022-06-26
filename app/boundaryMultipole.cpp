#include <iostream>
#include <fstream>
#include <armadillo>
#include "OData.h"
#include "BBH2D.h"
#include "BBH3D.h"
#include "MultiThread.h"
#include <chrono>

using namespace std;
using namespace arma;

void scanPolBBH2D(int nPoints, int dir, int * l, char* argv0, string name){
  OData out(argv0, name);
  double delta = 2/(double)nPoints;
  vector<vector<double>> data;
  vector<double> t1;
  vector<double> pol;
  BBH2D bbh2D(0,1);
  for(int i = 0; i <= nPoints; i++){
    if(i % (nPoints/10) == 0){
      cout << i << " %\r";
      cout.flush();
    }
    try{
      bbh2D.setIntraHop(i*delta);
      pol.push_back(bbh2D.getBoundPolarization(l,dir));
      t1.push_back(i*delta);
    }
    catch(const runtime_error & error){
      cout << "Singular matrix found at intracell hop = "  << i*delta << endl;
    }
  }
  cout << "100 %" << endl;
  data.push_back(t1);
  data.push_back(pol);
  out.data(data);
}

void scanQuadBBH3D(int nPoints, int dir, int * l, char* argv0, string name){
  OData out(argv0, name);
  double delta = 2/(double)nPoints;
  vector<vector<double>> data;
  vector<double> t1;
  vector<double> quad;
  BBH3D bbh3D(0,1);
  for(int i = 0; i <= nPoints; i++){
    if(i % (nPoints/10) == 0){
      cout << i << " %\r";
      cout.flush();
    }
    try{
      bbh3D.setIntraHop(i*delta);
      quad.push_back(bbh3D.getBoundQuadrupole(l,dir));
      t1.push_back(i*delta);
    }
    catch(const runtime_error & error){
      cout << "Singular matrix found at intracell hop = "  << i*delta << endl;
    }
  }
  cout << "100 %" << endl;
  data.push_back(t1);
  data.push_back(quad);
  out.data(data);
}

void threadPolBBH2D(int dir, int * l, vector<double> & res, vector<double> params){
  BBH2D bbh2D(0,1);
  double pol;
  for(int i = 0; i < params.size(); i++){
    try{
      bbh2D.setIntraHop(params[i]);
      pol = bbh2D.getBoundPolarization(l,dir);
      res.push_back(params[i]);
      res.push_back(pol);
    }
    catch(const runtime_error & error){
      cout << "Singular matrix found at intracell hop = "  << params[i] << endl;
    }
  }
}

void bbh2x1(vector<double> & res, vector<double> params){
  int l[2] = {50,100};
  threadPolBBH2D(0, l, res, params);
}

int main (int arc, char ** argv) {
  cout << "\n2D BBH" << endl;
  int l[2] = {50,100};
  /*
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
  MultiThread r1(bbh2x1, paramList, 8);
  r1.setFile(argv[0], "BoundaryPxBBH2D_50x100.dat");
  r1.run();

  cout << "\n3D BBH" << endl;
  int l2[3] = {10,10,10};
  BBH3D bbh(0.7,1);
  //cout << bbh.getBoundQuadrupole(l2, 0) << endl;
  //scanQuadBBH3D(100, 2, l2, argv[0], "BoundaryQxyBBH3D_10x10x10.dat");
  //scanQuadBBH3D(100, 1, l2, argv[0], "BoundaryQxzBBH3D_10x10x10.dat");
  //scanQuadBBH3D(100, 0, l2, argv[0], "BoundaryQyzBBH3D_10x10x10.dat");
  int l3[3] = {10,10,10};
  //scanQuadBBH3D(100, 2, l3, argv[0], "BoundaryQxyBBH3D_15x15x15.dat");
  //scanQuadBBH3D(100, 1, l3, argv[0], "BoundaryQxzBBH3D_15x15x15.dat");
  //scanQuadBBH3D(100, 0, l3, argv[0], "BoundaryQyzBBH3D_15x15x15.dat");

  return 0;
}
