#include <iostream>
#include "DisorderedSOTAI.h"
#include "ParallelMPI.h"

int sampPerJob = 40;
double m = 1.1;
int l[2] = {50,50};
int nStates = 50;

void lsrConstL(double * res, double * params){
  DisorderedSOTAI sotai(m);
  sotai.setSize(l);
  sotai.setW(params[0]);

  for(int i = 0; i < sampPerJob; i++){
    sotai.generateDisorder();
    try{
      res[i] = sotai.getLSR(nStates);
    }
    catch(const runtime_error & error){
      cout << "Matrix diagonalization failed for w = " << params[0] << endl;
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
  p.setFile(argv[0], "lsrSOTAI_L" + to_string(l[0]) + "_n" + to_string(nStates) + "_m1.1.dat");
  p.setJob(lsrConstL, sampPerJob);

  p.run();

}
