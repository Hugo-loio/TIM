#include <iostream>
#include "DisorderedSOTAI.h"
#include "ParallelMPI.h"

int sampPerJob = 1;
double m = 1.1;
int nStSamp = 5;

void loc(double * res, double * params){
  DisorderedSOTAI sotai(m);
  int l[2] = {(int)params[0], (int)params[0]};
  sotai.setSize(l);
  sotai.setW(params[1]);
  int nStates;
  double en = 0;
  if(params[1] == 0){
    en = 1e-5;
  }

  vector<double> resTemp;

  for(int i = 0; i < sampPerJob; i++){
    sotai.generateDisorder();
    try{
      nStates = 10*nStSamp;
      for(int e = 0; e < nStSamp; e++){
	resTemp.push_back(sotai.getIPR(nStates, en));
	resTemp.push_back(sotai.getLSR(nStates, en));
	nStates -= 10;
      }
      resTemp.push_back(sotai.getEnGap(en));
    }
    catch(const runtime_error & error){
      cout << "Matrix diagonalization failed for w = " << params[1] << " and L = " << params[0] << " at iteration " << i << endl;
      i--;
    }
  }

  for(int i = 0; i < resTemp.size(); i++){
    res[i] = resTemp[i];
  }
}


int main (int argc, char ** argv) {
  int sampMult = 200;
  int version = 0;
  if(argc > 1){
    version = stoi(argv[1]);
  }

  vector<vector<double>> paramList;
  int nPointsW = 100;
  int nPointsL = 20;
  for(int i = 2; i <= nPointsL; i++){
    for(int e = 0; e <= nPointsW; e++){
      if(i != 8){
	vector<double> param; 
	param.push_back(40 + 20*i);
	param.push_back(9*(double)e/(double)nPointsW);
	paramList.push_back(param);
      }
    }
  }

  ParallelMPI p(&argc, &argv);
  p.setSamples(sampMult);
  p.setParamList(paramList);
  p.setFile(argv[0], "locSOTAI_m1.1", version);
  p.setJob(loc, (1+2*nStSamp)*sampPerJob);
  p.setPrintEachSamp(false);
  p.run();
}
