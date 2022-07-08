#include <iostream>
#include <iomanip>
#include <thread>
#include "DisorderedSOTAI.h"
#include "OData.h"
#include "MultiThread.h"

void threadQuad(int * l, double m, int nSamples, vector<double> & res, vector <double> params){
  DisorderedSOTAI sotai(m);
  sotai.setSize(l);
  sotai.setW(params[0]);
  res.push_back(params[0]);

  for(int i = 0; i < nSamples; i++){
    sotai.generateDisorder();
    res.push_back(sotai.getQuadrupoleManyBody());
  }
}

void threadPol(int * l, double m, int nSamples, vector<double> & res, vector <double> params){
  DisorderedSOTAI sotai(m);
  sotai.setSize(l);
  sotai.setW(params[0]);
  res.push_back(params[0]);

  double polx;
  double poly;
  for(int i = 0; i < nSamples; i++){
    try{
      sotai.generateDisorder();
      polx = sotai.getBoundPolarization(0);
      poly = sotai.getBoundPolarization(1);
      res.push_back(polx);
      res.push_back(poly);
    }
    catch(const runtime_error & error){
      cout << "Singular matrix found at disorder weight = "  << params[0] << endl;
    }
  }
}

void run(void (*job) (vector<double> &, vector<double>), vector<vector<double>> & paramList, int threadNumber, char * argv0, string fileName, int version = 0){
  if(version != 0){
    fileName +=  "(" + to_string(version) + ").dat";
  }
  MultiThread r(job, paramList, threadNumber);
  r.setFile(argv0, fileName);
  r.run();
}

void quad1(vector<double> & res, vector<double> params){
  int l[2] = {10,10};
  threadQuad(l, 1.1, 40, res, params);
}

void quad2(vector<double> & res, vector<double> params){
  int l[2] = {20,20};
  threadQuad(l, 1.1, 40, res, params);
}

void quad3(vector<double> & res, vector<double> params){
  int l[2] = {40,40};
  threadQuad(l, 1.1, 40, res, params);
}

void pol1(vector<double> & res, vector<double> params){
  int l[2] = {10,10};
  threadPol(l, 1.1, 40, res, params);
}

void pol2(vector<double> & res, vector<double> params){
  int l[2] = {50,50};
  threadPol(l, 1.1, 40, res, params);
}

void pol3(vector<double> & res, vector<double> params){
  int l[2] = {100,100};
  threadPol(l, 1.1, 40, res, params);
}

void pol4(vector<double> & res, vector<double> params){
  int l[2] = {200,200};
  threadPol(l, 1.1, 40, res, params);
}

int main (int argc, char ** argv) {
  //SOTAI

  int threadNumber = 8;
  int version = 0;
  if(argc > 1){
    threadNumber = stoi(argv[1]);
  }
  if(argc > 2){
    version = stoi(argv[2]);
  }

  vector<vector<double>> paramList;
  int nPoints = 100;
  for(int i = 0; i <= nPoints; i++){
    vector<double> param; 
    param.push_back(9*(double)i/(double)nPoints);
    paramList.push_back(param);
  }

  //run(quad1, paramList, threadNumber, argv[0], "phaseDiagramSOTAI_10x10_m1.1",version);
  //run(quad2, paramList, threadNumber, argv[0], "phaseDiagramSOTAI_20x20_m1.1",version);
  //run(quad3, paramList, threadNumber, argv[0], "phaseDiagramSOTAI_40x40_m1.1",version);

  run(pol1, paramList, threadNumber, argv[0], "phaseDiagramSOTAIpol_10x10_m1.1",version);
  //run(pol2, paramList, threadNumber, argv[0], "phaseDiagramSOTAIpol_50x50_m1.1",version);
  //run(pol3, paramList, threadNumber, argv[0], "phaseDiagramSOTAIpol_100x100_m1.1",version);
  //run(pol4, paramList, threadNumber, argv[0], "phaseDiagramSOTAIpol_200x200_m1.1",version);
}
