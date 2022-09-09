#include <iostream>
#include "DisorderedSOTAI.h"
#include "ParallelMPI.h"

//int sampPerJob = 40;
double m = 1.1;
int l[2] = {100,100};
double w = 3.2;
int nMoments = 100;
int nRandVecs = 5;
double eMax = 10;

void dosConstW(double * res, double * params){
  DisorderedSOTAI sotai(m);
  sotai.setSize(l);
  sotai.setW(w);
  sotai.generateDisorder();

  for(int i = 0; i <= (int)params[0]; i++){
    res[2*i] = params[i+1];
    res[2*i+1] = sotai.getDOS(params[i+1], nMoments, nRandVecs, eMax);
  }
}

int main (int argc, char ** argv) {
  int sampMult = 40;

  vector<vector<double>> paramList1;
  int nPoints = 100;
  vector<double> param;
  param.push_back(nPoints);
  double deltaE = (2*eMax)/(double)nPoints;
  for(int i = 1; i <= nPoints + 1; i++){
    param.push_back(-eMax + deltaE*(i-1));
  }
  paramList1.push_back(param);

  ParallelMPI p(&argc, &argv);
  p.setSamples(sampMult);
  p.setParamList(paramList1);
  p.setFile(argv[0], "dosSOTAI_L" + to_string(l[0]) + "_nMu" + to_string(nMoments) + "_nR" + to_string(nRandVecs) + "_m1.1.dat");
  p.setJob(dosConstW, 2*(nPoints + 1));

  p.run();
}
