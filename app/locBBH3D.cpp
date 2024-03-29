#include <iostream>
#include "DisorderedBBH3D.h"
#include "ParallelMPI.h"

int sampPerJob = 1;
double m = 1.1;
int nStSamp = 5;

void loc(double * res, double * params){
  DisorderedBBH3D bbh3D(m);
  int l[3] = {(int)params[0], (int)params[0], (int)params[0]};
  bbh3D.setSize(l);
  bbh3D.setW(params[1]);
  int nStates;
  double en = 0;
  if(params[1] == 0){
    en = -0.01;
  }

  vector<double> resTemp;

  for(int i = 0; i < sampPerJob; i++){
    bbh3D.generateDisorder();
    try{
      nStates = 10*nStSamp;
      for(int e = 0; e < nStSamp; e++){
	resTemp.push_back(bbh3D.getIPR(nStates, en));
	resTemp.push_back(bbh3D.getLSR(nStates, en));
	nStates -= 10;
      }
      resTemp.push_back(bbh3D.getEnGap(en));
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
  for(int i = 9; i <= nPointsL; i++){
    for(int e = 1; e <= nPointsW; e++){
      vector<double> param; 
      param.push_back(4 + 2*i);
      param.push_back(9*(double)e/(double)nPointsW);
      paramList.push_back(param);
    }
  }

  ParallelMPI p(&argc, &argv);
  p.setSamples(sampMult);
  p.setParamList(paramList);
  p.setFile(argv[0], "locBBH3D_m1.1", version);
  p.setJob(loc, (1+2*nStSamp)*sampPerJob);
  p.setPrintEachSamp(false);
  p.run();
}
