#include <iostream>
#include "DisorderedBBH3D.h"
#include "ParallelMPI.h"
#include "AuxFunctions.h"

//int sampPerJob = 40;
double intra = 1.1;
int l[3] = {80,80,80};
double w = 2;
int nMoments = 4096;
int nRandVecs = 1;
int nPoints = 2000;

void dosConstW(double * res, double * params){
  DisorderedBBH3D bbh3d(intra);
  bbh3d.setSize(l);
  bbh3d.setW(w);
  bbh3d.generateDisorder();
  double eMax = w + 5;
  double deltaE = 2*eMax/(double)nPoints;

  for(int i = 0; i <= nPoints; i++){
    res[2*i] = -eMax + deltaE*i;
    res[2*i + 1] = bbh3d.getDOS(res[2*i], nMoments, nRandVecs, eMax);
  }
}

void dosE0(double * res, double * params){
  DisorderedBBH3D bbh3d(intra);
  bbh3d.setSize(l);
  bbh3d.setW(params[0]);
  bbh3d.generateDisorder();
  double eMax = (double)params[0] + 5;
  res[0] = bbh3d.getDOS(0, nMoments, nRandVecs, eMax);
}

int main (int argc, char ** argv) {
  if(argc > 1){
    w = stod(argv[1]);
  }
  int sampMult = 1;

  vector<vector<double>> paramList1;
  vector<double> param;
  for(int i = 0; i < sampMult; i++){
    paramList1.push_back(param);
  }

  vector<vector<double>> paramList2;
  int nPoints2 = 100;
  for(int i = 0; i <= nPoints2; i++){
    vector<double> param2;
    param2.push_back(9*(double)i/(double)nPoints2);
    paramList2.push_back(param2);
  }

  ParallelMPI p(&argc, &argv);
  if(argc > 1){
    p.setSamples(1);
    p.setParamList(paramList1);
    p.setFile(argv[0], "dosBBH3D_L" + to_string(l[0]) + "_w" + rmTrailZeros(to_string(w)) + "_nMu" + to_string(nMoments) + "_nR" + to_string(nRandVecs) + "_intra1.1");
    p.setJob(dosConstW, 2*(nPoints + 1));
    p.run();
  }
  else{
    p.setSamples(sampMult);
    p.setParamList(paramList2);
    p.setFile(argv[0], "dosBBH3D_L" + to_string(l[0]) + "_E0_nMu" + to_string(nMoments) + "_nR" + to_string(nRandVecs) + "_m1.1");
    p.setJob(dosE0, 1);
    p.run();
  }
}
