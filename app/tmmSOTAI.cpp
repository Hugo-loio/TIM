#include <iostream>
#include "DisorderedSOTAI.h"
#include "ParallelMPI.h"
#include "AuxFunctions.h"

double m = 1.1;
double w = 3.2;
double en = 0;
int qrIt = 10;
int l = 40;

void tmmConstW(double * res, double * params){
  DisorderedSOTAI sotai(m);
  sotai.setW(w);
  vector<double> tmm = sotai.getTMM(qrIt, en, params[0]);
  res[0] = tmm[0];
  res[1] = tmm[1];
}

void tmmConstL(double * res, double * params){
  DisorderedSOTAI sotai(m);
  sotai.setW(params[0]);
  vector<double> tmm = sotai.getTMM(qrIt, en, l);
  res[0] = tmm[0];
  res[1] = tmm[1];
}

int main (int argc, char ** argv) {
  if(argc > 1){
    en = stod(argv[1]);
  }

  vector<vector<double>> paramList1;
  int nPoints = 20;
  for(int i = 0; i <= nPoints; i++){
    vector<double> param; 
    param.push_back(12 + 4*i);
    paramList1.push_back(param);
  }

  vector<vector<double>> paramList2;
  for(int i = 59; i <= nPoints; i++){
    vector<double> param; 
    param.push_back(9*(double)i/(double)nPoints);
    paramList2.push_back(param);
  }

  ParallelMPI p(&argc, &argv);
  p.setParamList(paramList1);
  p.setFile(argv[0], "tmmSOTAI_E" + rmTrailZeros(to_string(en)) + "_w" + rmTrailZeros(to_string(w)) + "_m1.1");
  p.setJob(tmmConstW, 2);
  /*
     p.setParamList(paramList2);
     p.setFile(argv[0], "tmmSOTAI_E" + rmTrailZeros(to_string(en)) + "_L" + to_string(l) + "_m1.1");
     p.setJob(tmmConstL, 2);
     */
  p.run();

}
