#include <iostream>
#include "DisorderedSOTAI.h"
#include "ParallelMPI.h"
#include "AuxFunctions.h"

//int sampPerJob = 40;
double m = 1.1;
int l[2] = {100,100};
double w = 2;
int nMoments = 8192;
int nRandVecs = 1;
int nPoints = 2000;

void dos(double * res, double * params){
  DisorderedSOTAI sotai(m);
  sotai.setSize(l);
  sotai.setW(params[0]);
  sotai.generateDisorder();

  double eMax = (7/9)*params[0] + 4;

  for(int i = 0; i <= nPoints; i++){
    res[i] = sotai.getDOS(params[i+1], nMoments, nRandVecs, eMax);
  }
}

int main (int argc, char ** argv) {
  int sampMult = 200;

  /*
  vector<vector<double>> paramList;
  vector<double> param;
  double deltaE = (20)/(double)nPoints;
  for(int i = 0; i <= nPoints; i++){
    param.push_back(-10 + deltaE*i);
  }
  for(int i = 0; i < sampMult; i++){
    paramList1.push_back(param);
  }

  vector<vector<double>> paramList2;
  int nPoints2 = 100;
  for(int i = 10; i <= nPoints2; i++){
    vector<double> param2;
    param2.push_back(9*(double)i/(double)nPoints2);
    paramList2.push_back(param2);
  }

  double wVec[6] = {2.4, 2.8, 3.2, 3.6, 4, 9};

  ParallelMPI p(&argc, &argv);
  for(int i = 0; i < 0; i++){
    w = wVec[i];
    p.setSamples(1);
    p.setParamList(paramList1);
    p.setFile(argv[0], "dosSOTAI_L" + to_string(l[0]) + "_w" + rmTrailZeros(to_string(w)) + "_nMu" + to_string(nMoments) + "_nR" + to_string(nRandVecs) + "_m1.1");
    p.setJob(dosConstW, nPoints + 1);
    //p.setPrintEachSamp(true);
    p.run();
  }

  p.setSamples(40);
  p.setParamList(paramList2);
  p.setFile(argv[0], "dosSOTAI_L" + to_string(l[0]) + "_E0_nMu" + to_string(nMoments) + "_nR" + to_string(nRandVecs) + "_m1.1");
  p.setJob(dosE0PlusGap, 2);
  p.run();
  */
}
