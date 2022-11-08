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
  vector<double> tmm = bbh3d.getTMM(qrIt, params[0], params[1], dir);
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
  bool doCross = false;
  int part = 0;
  int nParts = 1;
  if(argc > 2){
    if(stoi(argv[1]) == 0){
      w = stod(argv[2]);
      if(argc > 3){
	dir = stod(argv[3]);
      }
    }
    else if(stoi(argv[1]) == 1){
      l = stoi(argv[2]);
      doConstW = false;
      if(argc > 3){
	dir = stod(argv[3]);
      }
    }
    else if(stoi(argv[1]) == 2){
      l = stoi(argv[2]);
      doConstW = false;
      doCross = true;
      if(argc > 3){
	dir = stod(argv[3]);
      }
    }
    if(argc > 5){
      nParts = stod(argv[4]);
      part = stod(argv[5]);
    }
  }

  vector<vector<double>> paramList1;
  int nPoints = 6;
  for(int e = 0; e <= 6; e++){
    for(int i = 5; i <= nPoints; i++){
      vector<double> param; 
      param.push_back(-3 + 0.5*e);
      param.push_back(2 + 2*i);
      paramList1.push_back(param);
    }
  }

  int nPoints2 = 150;
  vector<vector<double>> paramList2;
  for(int i = 0; i <= nPoints2; i++){
    vector<double> param; 
    param.push_back(9*(double)i/(double)100);
    paramList2.push_back(param);
  }

  int nPoints3 = 47;
  vector<vector<double>> paramList3;
  double wMin = 3.6;
  double wMax = 4.2;
  for(int i = 0; i <= nPoints3; i++){
    vector<double> param; 
    param.push_back(wMin + (wMax - wMin)*(double)i/(double)nPoints3);
    paramList3.push_back(param);
  }

  //Corrections for missing data
  if(doCross){
    if(l == 12 && dir == 2){
      paramList3.erase(paramList3.begin() , paramList3.begin() + 20);
      //printVec(paramList3);
    }
  }

  vector<vector<vector<double>>> paramLists;
  paramLists.push_back(paramList1);
  paramLists.push_back(paramList2);
  paramLists.push_back(paramList3);

  if(part != 0){
    for(int i = 0; i < paramLists.size(); i++){
      int n = paramLists[i].size();
      int begin = ((double)(part-1)/(double)nParts)*n;
      int end = ((double)part/(double)nParts)*n;
      paramLists[i] = vector<vector<double>>(&paramLists[i][begin], &paramLists[i][end]);
    }
  }

  /*
  ParallelMPI p(&argc, &argv);
  if(doConstW){
    p.setParamList(paramList1);
    p.setFile(argv[0], "tmmBBH3D_w" + rmTrailZeros(to_string(w)) + "_d" + to_string(dir) + "_m1.1", 0, part);
    p.setJob(tmmConstW, 2);
  }
  else{
    if(doCross){
      p.setParamList(paramList3);
      p.setFile(argv[0], "tmmBBH3D_E" + rmTrailZeros(to_string(en)) + "_L" + to_string(l) + "_d" + to_string(dir) + "_m1.1_cross", 0, part);
    }
    else{
      p.setParamList(paramList2);
      p.setFile(argv[0], "tmmBBH3D_E" + rmTrailZeros(to_string(en)) + "_L" + to_string(l) + "_d" + to_string(dir) + "_m1.1", 0, part);
    }
    p.setJob(tmmConstL, 2);
  }
  p.run();
  */
}
