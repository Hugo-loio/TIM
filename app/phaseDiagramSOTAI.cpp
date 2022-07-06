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
  if(argc == 2){
    threadNumber = stoi(argv[1]);
  }

  vector<vector<double>> paramList;
  int nPoints = 100;
  for(int i = 0; i <= nPoints; i++){
    vector<double> param; 
    param.push_back(9*(double)i/(double)nPoints);
    paramList.push_back(param);
  }

  /*
     MultiThread r1(quad1, paramList, threadNumber);
     r1.setFile(argv[0], "phaseDiagramSOTAI_10x10_m1.1.dat");
     r1.run();

     MultiThread r3(pol1, paramList, threadNumber);
     r3.setFile(argv[0], "phaseDiagramSOTAIpol_10x10_m1.1.dat");
     r3.run();

     MultiThread r2(quad2, paramList, threadNumber);
     r2.setFile(argv[0], "phaseDiagramSOTAI_20x20_m1.1.dat");
     r2.run();

     MultiThread r5(pol3, paramList, threadNumber);
     r5.setFile(argv[0], "phaseDiagramSOTAIpol_100x100_m1.1.dat");
     r5.run();

     MultiThread r7(quad3, paramList, threadNumber);
     r7.setFile(argv[0], "phaseDiagramSOTAI_40x40_m1.1.dat");
     r7.run();

     MultiThread r6(pol4, paramList, threadNumber);
     r6.setFile(argv[0], "phaseDiagramSOTAIpol_200x200_m1.1.dat");
     r6.run();
     */

     MultiThread r4(pol2, paramList, threadNumber);
     r4.setFile(argv[0], "phaseDiagramSOTAIpol_50x50_m1.1.dat");
     r4.run();

}
