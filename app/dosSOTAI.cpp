#include <iostream>
#include "DisorderedSOTAI.h"
#include "ParallelMPI.h"
#include "AuxFunctions.h"

//int sampPerJob = 40;
double m = 1.1;
int l[2] = {100,100};
double w = 3.2;
int nMoments = 100;
int nRandVecs = 1;
double eMax = 10;
int nPoints = 100;

void dosConstW(double * res, double * params){
  DisorderedSOTAI sotai(m);
  sotai.setSize(l);
  sotai.setW(w);
  sotai.generateDisorder();

  for(int i = 0; i <= nPoints; i++){
    res[i] = sotai.getDOS(params[i], nMoments, nRandVecs, eMax);
  }
}

int main (int argc, char ** argv) {
  int sampMult = 40;

  vector<vector<double>> paramList1;
  vector<double> param;
  double deltaE = (18)/(double)nPoints;
  for(int i = 0; i <= nPoints; i++){
    param.push_back(-9 + deltaE*i);
  }
  for(int i = 0; i < sampMult; i++){
    paramList1.push_back(param);
  }

  ParallelMPI p(&argc, &argv);
  p.setSamples(1);
  p.setParamList(paramList1);
  p.setFile(argv[0], "dosSOTAI_L" + to_string(l[0]) + "_w" + rmTrailZeros(to_string(w)) + "_nMu" + to_string(nMoments) + "_nR" + to_string(nRandVecs) + "_m1.1.dat");
  p.setJob(dosConstW, nPoints + 1);
  p.setPrintEachSamp(true);

  p.run();
}
