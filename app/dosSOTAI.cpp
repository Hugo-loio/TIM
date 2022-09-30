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

void dosConstW(double * res, double * params){
  DisorderedSOTAI sotai(m);
  sotai.setSize(l);
  sotai.setW(w);
  sotai.generateDisorder();
  double eMax = w*7/9 + 4;
  double deltaE = 2*eMax/(double)nPoints;

  for(int i = 0; i <= nPoints; i++){
    res[2*i] = -eMax + deltaE*i;
    res[2*i + 1] = sotai.getDOS(res[2*i], nMoments, nRandVecs, eMax);
  }
}

void dosE0(double * res, double * params){
  DisorderedSOTAI sotai(m);
  sotai.setSize(l);
  sotai.setW(params[0]);
  sotai.generateDisorder();
  double eMax = (double)params[0]*7/9 + 4;
  res[0] = sotai.getDOS(0, nMoments, nRandVecs, eMax);
}

int main (int argc, char ** argv) {
  int sampMult = 200;

  vector<vector<double>> paramList1;
  vector<double> param;
  for(int i = 0; i < sampMult; i++){
    paramList1.push_back(param);
  }

  vector<vector<double>> paramList2;
  int nPoints2 = 100;
  for(int i = 75; i <= nPoints2; i++){
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
    p.setJob(dosConstW, 2*(nPoints + 1));
    p.run();
  }

  p.setSamples(sampMult);
  p.setParamList(paramList2);
  p.setFile(argv[0], "dosSOTAI_L" + to_string(l[0]) + "_E0_nMu" + to_string(nMoments) + "_nR" + to_string(nRandVecs) + "_m1.1");
  p.setJob(dosE0, 1);
  p.run();
}
