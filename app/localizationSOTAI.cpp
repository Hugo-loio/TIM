#include <iostream>
#include <iomanip>
#include <thread>
#include "DisorderedSOTAI.h"
#include "OData.h"
#include "MultiThread.h"
#include "AuxFunctions.h"

void ipr(int * l, double m, int nSamples, int nStates, vector<double> & res, vector <double> params){
  DisorderedSOTAI sotai(m);
  sotai.setSize(l);
  sotai.setW(params[0]);
  res.push_back(params[0]);

  for(int i = 0; i < nSamples; i++){
    sotai.generateDisorder();
    res.push_back(sotai.getIPR(nStates));
  }
}

void ipr1(vector<double> & res, vector<double> params){
  int l[2] = {20,20};
  ipr(l, 1.1, 5, 10, res, params);
}

void ipr2(vector<double> & res, vector<double> params){
  int l[2] = {50,50};
  ipr(l, 1.1, 5, 10, res, params);
}

void ipr3(vector<double> & res, vector<double> params){
  int l[2] = {200,200};
  ipr(l, 1.1, 5, 10, res, params);
}

int main (int argc, char ** argv) {
  int threadNumber = 8;
  int version = 0;
  int part = 0;
  if(argc > 1){
    threadNumber = stoi(argv[1]);
  }
  if(argc > 2){
    version = stoi(argv[2]);
  }
  if(argc > 3){
    part = stoi(argv[3]);
  }

  vector<vector<double>> paramList;
  int nPoints = 100;
  for(int i = 0; i <= nPoints; i++){
    vector<double> param; 
    param.push_back(9*(double)i/(double)nPoints);
    paramList.push_back(param);
  }

  int nParts = 4;
  if(part != 0){
    if(part != 1){
      paramList.erase(paramList.begin(), paramList.begin() + (nPoints/nParts)*(part -1));
    }
    if(part != nParts){
      paramList.erase(paramList.begin() + nPoints/nParts, paramList.end());
    }
  }

  //runSingleThread(ipr1, paramList, argv[0], "iprSOTAI_L20_n10_m1.1", version, part);
  //runSingleThread(ipr2, paramList, argv[0], "iprSOTAI_L50_n10_m1.1", version, part);
  runSingleThread(ipr3, paramList, argv[0], "iprSOTAI_L200_n10_m1.1", version, part);
}
