#include <iostream>
#include "DisorderedBBH3D.h"
#include "ParallelMPI.h"
#include "AuxFunctions.h"
#include <cmath>

//int sampPerJob = 40;
double m = 1.1;
int l[3] = {10,10,10};
double w = 3;
int nMoments = 1024;
int nPoints = 10;
double en = 0;

void ldosEdge(double * res, double * params){
  DisorderedBBH3D bbh3d(m);
  bbh3d.setSize(l);
  bbh3d.setW(w);
  bbh3d.generateDisorder();
  double eMax = w + 5;
  int count = 0;
  int n[3] = {0,0,0};

  for(int i = 0; i < l[1] && i < nPoints; i++){
    n[0] = i;
    //cout << e << " " << i << " " << j << endl;
    //cout << count << endl;
    res[count++] = n[0];
    res[count++] = bbh3d.getLDOS(n, en, nMoments, eMax);
  }
}

void ldosDiag(double * res, double * params){
  DisorderedBBH3D bbh3d(m);
  bbh3d.setSize(l);
  bbh3d.setW(w);
  bbh3d.generateDisorder();
  double eMax = w + 5;
  int count = 0;
  int n[3] = {0,0,0};

  for(int i = 0; i < l[1] && i < nPoints; i++){
    n[0] = i;
    n[1] = i;
    n[2] = i;
    res[count++] = sqrt(3*i*i);
    res[count++] = bbh3d.getLDOS(n, en, nMoments, eMax);
  }
}

int main (int argc, char ** argv) {
  int sampMult = 200;
  int mode = 0;
  if(argc > 1){
    w = stod(argv[1]);
    if(argc > 2){
      l[0] = stoi(argv[2]);
      l[1] = l[0];
      l[2] = l[0];
      if(argc > 3){
	nMoments = stoi(argv[3]);
	if (argc > 4){
	  mode = stoi(argv[4]);
	}
      }
    }
  }

  vector<vector<double>> paramList1;
  vector<double> param;
  for(int i = 0; i < sampMult; i++){
    paramList1.push_back(param);
  }

  int resSize = 2;
  if(l[0] < nPoints){
    resSize *= l[0];
  }
  else{
    resSize *= nPoints;
  }

  ParallelMPI p(&argc, &argv);
  //p.setSamples(sampMult);
  p.setParamList(paramList1);
  if(mode == 0){
    p.setFile(argv[0], "ldosBBH3D_L" + to_string(l[0]) + "_w" + rmTrailZeros(to_string(w)) + "_E0_nMu" + to_string(nMoments) + "_m1.1_edge");
    p.setJob(ldosEdge, resSize);
  }
  else if(mode == 1){
    p.setFile(argv[0], "ldosBBH3D_L" + to_string(l[0]) + "_w" + rmTrailZeros(to_string(w)) + "_E0_nMu" + to_string(nMoments) + "_m1.1_diag");
    p.setJob(ldosDiag, resSize);
  }
  p.run();
}
