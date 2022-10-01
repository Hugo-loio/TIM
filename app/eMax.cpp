#include <iostream>
#include "DisorderedSOTAI.h"
#include "DisorderedBBH3D.h"
#include "ParallelMPI.h"
#include "AuxFunctions.h"

int nHam = 50;
double m = 1.1;
int l = 6;

void eMaxSOTAI(double * res, double * params){
  DisorderedSOTAI sotai(m);
  int lVec[2] = {l,l};
  sotai.setSize(lVec);
  sotai.setW(params[0]);
  double max = 0;
  double maxTemp;

  for(int i = 0; i <= nHam; i++){
    sotai.generateDisorder();
    maxTemp = sotai.getMaxE();
    if(maxTemp > max){
      max = maxTemp;
    }
  }
  res[0] = max;
}

void eMaxBBH3D(double * res, double * params){
  DisorderedBBH3D bbh(m);
  int lVec[3] = {l,l,l};
  bbh.setSize(lVec);
  bbh.setW(params[0]);
  double max = 0;
  double maxTemp;

  for(int i = 0; i <= nHam; i++){
    bbh.generateDisorder();
    maxTemp = bbh.getMaxE();
    if(maxTemp > max){
      max = maxTemp;
    }
  }
  res[0] = max;
}

int main (int argc, char ** argv) {
  vector<vector<double>> paramList;
  int nPoints = 100;
  for(int i = 0; i <= nPoints; i++){
    vector<double> param;
    param.push_back(9*(double)i/(double)nPoints);
    paramList.push_back(param);
  }


  ParallelMPI p(&argc, &argv);
  p.setSamples(1);
  p.setParamList(paramList);
  //p.setFile(argv[0], "eMaxSOTAI_L" + to_string(l) + "_nHam" + to_string(nHam) + "_m1.1");
  //p.setJob(eMaxSOTAI, 1);
  p.setFile(argv[0], "eMaxBBH3D_L" + to_string(l) + "_nHam" + to_string(nHam) + "_m1.1");
  p.setJob(eMaxBBH3D, 1);
  p.run();
}
