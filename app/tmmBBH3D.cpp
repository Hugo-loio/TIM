#include <iostream>
#include "DisorderedBBH3D.h"
#include "ParallelMPI.h"
#include "AuxFunctions.h"

double m = 1.1;
double w = 3.2;
double en = 0;
int qrIt = 10;
int l = 40;
int dir = 2;

void tmmConstW(double * res, double * params){
  DisorderedBBH3D bbh3d(m);
  bbh3d.setW(w);
  vector<double> tmm = bbh3d.getTMM(qrIt, en, params[0], dir);
  res[0] = tmm[0];
  res[1] = tmm[1];
}

void tmmConstL(double * res, double * params){
  DisorderedBBH3D bbh3d(m);
  bbh3d.setW(params[0]);
  vector<double> tmm = bbh3d.getTMM(qrIt, en, l, dir);
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
      if(argc > 4){
	dir = stod(argv[4]);
      }
    }
    else if(stoi(argv[1]) == 1){
      l = stoi(argv[2]);
      doConstW = false;
      if(argc > 3){
	dir = stod(argv[3]);
      }
    }
  }

  vector<vector<double>> paramList1;
  int nPoints = 10;
  for(int i = 0; i <= nPoints; i++){
    vector<double> param; 
    param.push_back(2 + 2*i);
    paramList1.push_back(param);
  }

  int nPoints2 = 150;
  vector<vector<double>> paramList2;
  for(int i = 0; i <= nPoints2; i++){
    vector<double> param; 
    param.push_back(9*(double)i/(double)100);
    paramList2.push_back(param);
  }

  //Corrections for missing data
  if(doConstW){
  }
  else{
    if(l == 12 && dir == 2){
      paramList2.erase(paramList2.begin() + 55, paramList2.begin() + 57);
      paramList2.erase(paramList2.begin() + 30, paramList2.begin() + 54);
      paramList2.erase(paramList2.begin(), paramList2.begin() + 28);
      //printVec(paramList2);
    }
    if(l == 14 && dir == 2){
      paramList2.erase(paramList2.begin() + 47, paramList2.begin() + 51);
      paramList2.erase(paramList2.begin() + 31, paramList2.begin() + 46);
      //printVec(paramList2);
    }
  }


  ParallelMPI p(&argc, &argv);
  if(doConstW){
    p.setParamList(paramList1);
    p.setFile(argv[0], "tmmBBH3D_E" + rmTrailZeros(to_string(en)) + "_w" + rmTrailZeros(to_string(w)) + "_d" + to_string(dir) + "_m1.1");
    p.setJob(tmmConstW, 2);
  }
  else{
    p.setParamList(paramList2);
    p.setFile(argv[0], "tmmBBH3D_E" + rmTrailZeros(to_string(en)) + "_L" + to_string(l) + "_d" + to_string(dir) + "_m1.1");
    p.setJob(tmmConstL, 2);
  }
  p.run();
}
