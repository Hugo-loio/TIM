#include <iostream>
#include "DisorderedBBH3D.h"
#include "ParallelMPI.h"
#include "AuxFunctions.h"
#include <cmath>

//int sampPerJob = 40;
double m = 1.1;
int l[3] = {10,10,10};
double en = 0;
int nStates = 50;

void spectrum(double * res, double * params){
  DisorderedBBH3D bbh3d(m);
  bbh3d.setSize(l);
  bbh3d.setW(params[0]);
  bbh3d.generateDisorder();

  vec en = bbh3d.spectrumE0(nStates);
  for(int i = 0; i < size(en)[0]; i++){
    res[i] = en[i];
  }
}

int main (int argc, char ** argv) {
  int sampMult = 200;
  int mode = 0;
  if(argc > 1){
    l[0] = stoi(argv[1]);
    l[1] = l[0];
    l[2] = l[0];
    if(argc > 2){
      mode = stoi(argv[2]);
    }
  }

  vector<vector<double>> paramList;
  int nPointsW = 100;
  for(int e = 0; e <= nPointsW; e++){
    vector<double> param; 
    param.push_back(9*(double)e/(double)nPointsW);
    paramList.push_back(param);
  }

  vector<vector<double>> paramList2;
  double startW = 3.5;
  double endW = 4;
  for(int e = 0; e <= nPointsW; e++){
    vector<double> param; 
    param.push_back(startW+(double)e*(endW-startW)/(double)nPointsW);
    paramList2.push_back(param);
  }

  ParallelMPI p(&argc, &argv);
  p.setSamples(sampMult);
  if(mode == 0){
    p.setParamList(paramList);
    p.setFile(argv[0], "spectrumE0BBH3D_L" + to_string(l[0]) + "_m1.1");
  }
  else if(mode == 1){
    p.setParamList(paramList2);
    p.setFile(argv[0], "spectrumE0BBH3D_L" + to_string(l[0]) + "_m1.1_cross");
  }
  p.setJob(spectrum, nStates);
  p.run();
}
