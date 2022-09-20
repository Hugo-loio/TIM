#include <iostream>
#include "DisorderedSOTAI.h"
#include "ParallelMPI.h"
#include "AuxFunctions.h"

int sampPerJob = 200;
double m = 1.1;
int l[2] = {100,100};
int nStates = 10;

void gapConstL(double * res, double * params){
  DisorderedSOTAI sotai(m);
  sotai.setSize(l);
  sotai.setW(params[0]);
  double en = 0;
  if(params[0] == 0){
    en = 1e-5;
  }

  for(int i = 0; i < sampPerJob; i++){
    sotai.generateDisorder();
    try{
      res[i] = sotai.getEnGap(en);
    }
    catch(const runtime_error & error){
      cout << "Matrix diagonalization failed for w = " << params[0] << " at iteration " << i << endl;
      i--;
    }
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

  ParallelMPI p(&argc, &argv);
  p.setSamples(sampMult);
  p.setParamList(paramList1);
  p.setFile(argv[0], "enGapSOTAI_L" + to_string(l[0]) + "_m1.1");
  p.setJob(gapConstL, sampPerJob);
  p.setPrintEachSamp(true);
  p.run();
}
