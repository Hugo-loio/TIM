#include <iostream>
#include "Anderson3D.h"
#include "ParallelMPI.h"
#include "AuxFunctions.h"

double m = 1.1;
double w = 3.2;
double en = 0;
int qrIt = 10;
int l = 40;

void tmmConstW(double * res, double * params){
  Anderson3D and3d(m);
  and3d.setW(w);
  vector<double> tmm = and3d.getTMM(qrIt, en, params[0]);
  res[0] = tmm[0];
  res[1] = tmm[1];
}

void tmmConstL(double * res, double * params){
  Anderson3D and3d(m);
  and3d.setW(params[0]);
  vector<double> tmm = and3d.getTMM(qrIt, en, l);
  res[0] = tmm[0];
  res[1] = tmm[1];
}

void tmm(double * res, double * params){
  Anderson3D and3d(m);
  and3d.setW(params[0]);
  //cout << params[0] << " " << params[1] << endl;
  vector<double> tmm = and3d.getTMM(qrIt, en, params[1]);
  res[0] = tmm[0];
  res[1] = tmm[1];
}

int main (int argc, char ** argv) {
  bool doConstW = true;
  if(argc > 2){
    if(stoi(argv[1]) == 0){
      w = stod(argv[2]);
      if(argc > 3){
	en = stod(argv[3]);
      }
    }
    else if(stoi(argv[1]) == 1){
      l = stoi(argv[2]);
      doConstW = false;
    }
  }

  vector<vector<double>> paramList1;
  int nPoints = 8;
  for(int i = 0; i <= nPoints; i++){
    vector<double> param; 
    param.push_back(1 + i);
    paramList1.push_back(param);
  }

  int nPoints2 = 150;
  vector<vector<double>> paramList2;
  for(int i = 0; i <= nPoints2; i++){
    vector<double> param; 
    param.push_back(9*(double)i/(double)100);
    paramList2.push_back(param);
  }

  vector<vector<double>> paramList3;
  int nW = 20;
  for(int i = 0; i <= nW; i++){
    for(int e = 0; e <= nPoints; e++){
      vector<double> param; 
      param.push_back(10 + i);
      param.push_back(1 + e);
      paramList3.push_back(param);
    }
  }

  ParallelMPI p(&argc, &argv);
  p.setParamList(paramList3);
  p.setFile(argv[0], "tmmAnderson3D_E" + rmTrailZeros(to_string(en)) + "_t1");
  p.setJob(tmm, 2);
  /*
     if(doConstW){
     p.setParamList(paramList1);
     p.setFile(argv[0], "tmmAnderson3D_E" + rmTrailZeros(to_string(en)) + "_w" + rmTrailZeros(to_string(w)) + "_m1.1");
     p.setJob(tmmConstW, 2);
     }
     else{
     p.setParamList(paramList2);
     p.setFile(argv[0], "tmmAnderson3D_E" + rmTrailZeros(to_string(en)) + "_L" + to_string(l) + "_m1.1");
     p.setJob(tmmConstL, 2);
     }
     */
  p.run();
}
