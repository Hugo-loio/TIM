#include <iostream>
#include "DisorderedSOTAI.h"
#include "ParallelMPI.h"
#include "AuxFunctions.h"

//int sampPerJob = 40;
double m = 1.1;
int l[2] = {10,10};
double w = 3;
int nMoments = 2048;
double en = 0;

void ldos(double * res, double * params){
  DisorderedSOTAI sotai(m);
  sotai.setSize(l);
  sotai.setW(w);
  sotai.generateDisorder();
  double eMax = w*7/9 + 4;
  int count = 0;
  int n[2] = {0,0};

  for(int i = 0; i < l[1]; i++){
    n[1] = i;
    for(int e = 0; e < l[0]; e++){
      n[0] = e;
      res[count++] = n[0];
      res[count++] = n[1];
      res[count++] = sotai.getLDOS(n, en, nMoments, eMax);
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
  p.setFile(argv[0], "ldosSOTAI_L" + to_string(l[0]) + "_w" + rmTrailZeros(to_string(w)) + "_E0_nMu" + to_string(nMoments) + "_m1.1");
  p.setJob(ldos, 3*l[0]*l[1]);
  p.run();
}
