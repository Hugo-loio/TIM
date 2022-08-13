#include <iostream>
#include <iomanip>
#include <thread>
#include "DisorderedSOTAI.h"
#include "OData.h"
#include "MultiThread.h"
#include "AuxFunctions.h"

void threadQuad(int * l, double m, vector<double> & res, vector <double> params){
  DisorderedSOTAI sotai(m);
  sotai.setSize(l);
  sotai.setW(params[0]);
  sotai.generateDisorder();
  res.push_back(sotai.getQuadrupoleManyBody());
}

void threadPol(int * l, double m, vector<double> & res, vector <double> params){
  DisorderedSOTAI sotai(m);
  sotai.setSize(l);
  sotai.setW(params[0]);

  try{
    sotai.generateDisorder();
    double polx = sotai.getBoundPolarization(0);
    double poly = sotai.getBoundPolarization(1);
    res.push_back(polx);
    res.push_back(poly);
  }
  catch(const runtime_error & error){
    cout << "Singular matrix found at disorder weight = "  << params[0] << endl;
  }
}

void quad1(vector<double> & res, vector<double> params){
  int l[2] = {10,10};
  threadQuad(l, 1.1, res, params);
}

void quad2(vector<double> & res, vector<double> params){
  int l[2] = {20,20};
  threadQuad(l, 1.1, res, params);
}

void quad3(vector<double> & res, vector<double> params){
  int l[2] = {40,40};
  threadQuad(l, 1.1, res, params);
}

void quad4(vector<double> & res, vector<double> params){
  int l[2] = {60,60};
  threadQuad(l, 1.1, res, params);
}

void quad5(vector<double> & res, vector<double> params){
  int l[2] = {80,80};
  threadQuad(l, 1.1, res, params);
}

void pol1(vector<double> & res, vector<double> params){
  int l[2] = {10,10};
  threadPol(l, 1.1, res, params);
}

void pol2(vector<double> & res, vector<double> params){
  int l[2] = {50,50};
  threadPol(l, 1.1, res, params);
}

void pol3(vector<double> & res, vector<double> params){
  int l[2] = {100,100};
  threadPol(l, 1.1, res, params);
}

void pol4(vector<double> & res, vector<double> params){
  int l[2] = {200,200};
  threadPol(l, 1.1, res, params);
}

void pol5(vector<double> & res, vector<double> params){
  int l[2] = {300,300};
  threadPol(l, 1.1, res, params);
}

void pol6(vector<double> & res, vector<double> params){
  int l[2] = {400,400};
  threadPol(l, 1.1, res, params);
}

void pol7(vector<double> & res, vector<double> params){
  int l[2] = {500,500};
  threadPol(l, 1.1, res, params);
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

  //run(quad1, paramList, threadNumber, argv[0], "phaseDiagramSOTAI_10x10_m1.1", 40,version, part);
  //run(quad2, paramList, threadNumber, argv[0], "phaseDiagramSOTAI_20x20_m1.1", 40, version, part);
  //run(quad3, paramList, threadNumber, argv[0], "phaseDiagramSOTAI_40x40_m1.1", 40, version, part);
  //run(quad4, paramList, threadNumber, argv[0], "phaseDiagramSOTAI_60x60_m1.1", 40, version, part);
  //run(quad5, paramList, threadNumber, argv[0], "phaseDiagramSOTAI_80x80_m1.1", 40, version, part);

  run(pol1, paramList, threadNumber, argv[0], "phaseDiagramSOTAIpol_10x10_m1.1", 40, version, part);
  //run(pol2, paramList, threadNumber, argv[0], "phaseDiagramSOTAIpol_50x50_m1.1", 40, version, part);
  //run(pol3, paramList, threadNumber, argv[0], "phaseDiagramSOTAIpol_100x100_m1.1", 40, version, part);
  //run(pol4, paramList, threadNumber, argv[0], "phaseDiagramSOTAIpol_200x200_m1.1", 40, version, part);
  //run(pol5, paramList, threadNumber, argv[0], "phaseDiagramSOTAIpol_300x300_m1.1", 40, version, part);
  //run(pol6, paramList, threadNumber, argv[0], "phaseDiagramSOTAIpol_400x400_m1.1", 40, version, part);
  //run(pol7, paramList, threadNumber, argv[0], "phaseDiagramSOTAIpol_500x500_m1.1", 40, version, part);
}
