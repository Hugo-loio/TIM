#include <iostream>
#include "DisorderedSOTAI.h"
#include "ParallelMPI.h"

int sampPerJob = 10;
double m = 1.1;
int l[2] = {20, 20};
int nStates = 10;

void ipr(double * res, double * params){
  DisorderedSOTAI sotai(m);
  sotai.setSize(l);
  sotai.setW(params[0]);

  for(int i = 0; i < sampPerJob; i++){
    sotai.generateDisorder();
    res[i] = sotai.getIPR(nStates);
  }
}

int main (int argc, char ** argv) {
  int sampMult = 4;

  vector<vector<double>> paramList;
  int nPoints = 100;
  for(int i = 0; i <= nPoints; i++){
    vector<double> param; 
    param.push_back(9*(double)i/(double)nPoints);
    paramList.push_back(param);
  }

  ParallelMPI p(&argc, &argv);
  p.setSamples(sampMult);
  p.setParamList(paramList);
  p.setFile(argv[0], "iprSOTAI_L" + to_string(l[0]) + "_n" + to_string(nStates) + "_m1.1.dat");
  p.setJob(ipr, sampPerJob);
  p.run();

}
