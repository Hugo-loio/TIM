#include <iostream>
#include "DisorderedBBH3D.h"
#include "ParallelMPI.h"
#include "AuxFunctions.h"

//int sampPerJob = 40;
double m = 1.1;
int l[3] = {4,4,4};
double w = 3;
int nMoments = 1024;
double en = 0;

void ldos(double * res, double * params){
  DisorderedBBH3D bbh3d(m);
  bbh3d.setSize(l);
  bbh3d.setW(w);
  bbh3d.generateDisorder();
  double eMax = w + 5;
  int count = 0;
  int n[3] = {0,0,0};

  for(int j = 0; j < l[2] || j < 5; j++){
    n[2] = j;
    for(int i = 0; i < l[1] || j < 5; i++){
      n[1] = i;
      for(int e = 0; e < l[0] || j < 5; e++){
	n[0] = e;
	res[count++] = n[0];
	res[count++] = n[1];
	res[count++] = n[2];
	res[count++] = bbh3d.getLDOS(n, en, nMoments, eMax);
      }
    }
  }
}

int main (int argc, char ** argv) {
  int sampMult = 200;
  if(argc > 1){
    w = stod(argv[1]);
    if(argc > 2){
      l[0] = stoi(argv[2]);
      l[1] = l[0];
      l[2] = l[0];
      if(argc > 3){
	nMoments = stoi(argv[3]);
      }
    }
  }

  vector<vector<double>> paramList1;
  vector<double> param;
  for(int i = 0; i < sampMult; i++){
    paramList1.push_back(param);
  }

  ParallelMPI p(&argc, &argv);
  //p.setSamples(sampMult);
  p.setParamList(paramList1);
  p.setFile(argv[0], "ldosBBH3D_L" + to_string(l[0]) + "_w" + rmTrailZeros(to_string(w)) + "_E0_nMu" + to_string(nMoments) + "_m1.1");
  p.setJob(ldos, 4*l[0]*l[1]*l[2]);
  p.run();
}
