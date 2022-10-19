#include <iostream>
#include "DisorderedSOTAI.h"
#include "ParallelMPI.h"
#include "AuxFunctions.h"

int sampPerJob = 1;
double m = 1.1;
double w = 3.2;

void fractal(double * res, double * params){
  DisorderedSOTAI sotai(m);
  int l[2] = {(int)params[0], (int)params[0]};
  sotai.setSize(l);
  sotai.setW(w);

  for(int i = 0; i < sampPerJob; i++){
    sotai.generateDisorder();
    try{
      res[i] = sotai.getIPR(10, params[1]);
    }
    catch(const runtime_error & error){
      cout << "Matrix diagonalization failed for en = " << params[1] << " and L = " << params[0] << " at iteration " << i << endl;
      i--;
    }
  }
}


int main (int argc, char ** argv) {
  int sampMult = 200;
  int version = 0;
  if(argc > 1){
    w = stod(argv[1]);
    if(argc > 2){
      version = stoi(argv[2]);
    }
  }

  vector<vector<double>> paramList;
  int nPointsE = 20;
  int nPointsL = 20;
  for(int i = 0; i <= nPointsL; i++){
    if((160 + 4*i) % 20 != 0){
      for(int e = 0; e <= nPointsE; e++){
	vector<double> param; 
	param.push_back(160 + 4*i);
	param.push_back(-3 + 3*(double)e/(double)nPointsE);
	paramList.push_back(param);
      }
    }
  }

  ParallelMPI p(&argc, &argv);
  p.setSamples(sampMult);
  p.setParamList(paramList);
  p.setFile(argv[0], "fractalSOTAI_w" + rmTrailZeros(to_string(w)) + "_m1.1", version);
  p.setJob(fractal, sampPerJob);
  p.setPrintEachSamp(false);
  p.run();
}
