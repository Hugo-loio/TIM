#include <iostream>
#include "DisorderedSOTAI.h"
#include "ParallelMPI.h"

int sampPerJob = 40;
double m = 1.1;
int l[2] = {50,50};
int nStates = 10;

void iprConstL(double * res, double * params){
  DisorderedSOTAI sotai(m);
  sotai.setSize(l);
  sotai.setW(params[0]);

  for(int i = 0; i < sampPerJob; i++){
    sotai.generateDisorder();
    res[i] = sotai.getIPR(nStates);
  }
}

void ipr(double * res, double * params){
  DisorderedSOTAI sotai(m);
  int l2[2] = {(int)params[0], (int)params[0]};
  sotai.setSize(l2);
  sotai.setW(params[1]);

  for(int i = 0; i < sampPerJob; i++){
    sotai.generateDisorder();
    res[i] = sotai.getIPR(nStates);
  }
}

int main (int argc, char ** argv) {
  int sampMult = 1;

  vector<vector<double>> paramList1;
  int nPoints = 100;
  for(int i = 0; i <= nPoints; i++){
    vector<double> param; 
    param.push_back(9*(double)i/(double)nPoints);
    paramList1.push_back(param);
  }

  int nPointsFit = 15;
  vector<vector<double>> paramList2;
  for(int i = 0; i <= nPointsFit; i++){
    for(int e = 0; e <= nPoints; e++){
      vector<double> param;
      param.push_back(90 + 10*i);
      param.push_back(9*(double)e/(double)nPoints);
      paramList2.push_back(param);
    }
  }

  ParallelMPI p(&argc, &argv);
  p.setSamples(sampMult);
  /*
     p.setParamList(paramList1);
     p.setFile(argv[0], "iprSOTAI_L" + to_string(l[0]) + "_n" + to_string(nStates) + "_m1.1.dat");
     p.setJob(iprConstL, sampPerJob);
     */
  p.setParamList(paramList2);
  p.setFile(argv[0], "iprSOTAI_n" + to_string(nStates) + "_m1.1.dat");
  p.setJob(ipr, sampPerJob);

  p.run();

}
