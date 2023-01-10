#include <iostream>
#include "DisorderedBBH3D.h"
#include "ParallelMPI.h"
#include "AuxFunctions.h"
#include <cmath>

//int sampPerJob = 40;
double m = 1.1;
int l[3] = {10,10,10};
double w = 3;
double en = 0;

void probDens(double * res, double * params){
  DisorderedBBH3D bbh3d(m);
  bbh3d.setSize(l);
  bbh3d.setW(w);
  bbh3d.generateDisorder();
  double eMax = w + 5;
  int count = 0;
  int n[3] = {0,0,0};

  for(int i = 0; i < l[0]; i++){
    for(int j = 0; j < l[1]; j++){
      for(int k = 0; k < l[2]; k++){
	n[0] = i;
	n[1] = j;
	n[2] = k;
	res[count++] = n[0];
	res[count++] = n[1];
	res[count++] = n[2];
	res[count++] = bbh3d.probDensE0(n);
      }
    }
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
  p.setFile(argv[0], "propDensBBH3D_L" + to_string(l[0]) + "_w" + rmTrailZeros(to_string(w)) + "_E0_m1.1");
  p.setJob(probDens, 4*l[0]*l[1]*l[2]);
  p.run();
}
