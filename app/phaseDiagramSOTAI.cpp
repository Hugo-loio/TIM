#include <iostream>
#include <iomanip>
#include <thread>
#include "DisorderedSOTAI.h"
#include "OData.h"
#include "MultiThread.h"

//Multithreading solution
void threadCalcQuad(int nAvg, int nThreads, int threadIndex, int * l, OData & out){
  int nPoints = 21;
  double m = 1.1;
  int pointsPerThread = nPoints/nThreads;
  int addPoints = 0;
  if(nPoints % nThreads > threadIndex){
    pointsPerThread++;
  }
  else{
    addPoints = nPoints % nThreads;
  }
  int nStart = pointsPerThread*threadIndex + addPoints;
  int nEnd = pointsPerThread*(threadIndex + 1) + addPoints;
  double delta = 9/(double)(nPoints - 1);

  vector<vector<double>> data;
  vector<double> w;
  vector<vector<double>> quad;
  DisorderedSOTAI sotai(m);
  sotai.setSize(l);

  for(int e = 0; e < nAvg; e++){
    quad.push_back(w);
  }

  for(int i = nStart; i < nEnd; i++){
    w.push_back(i*delta);
    sotai.setW(i*delta);
    cout << i << " " << threadIndex << endl;
    for(int e = 0; e < nAvg; e++){
      sotai.generateDisorder();
      quad[e].push_back(sotai.getQuadrupoleManyBody());
    }
  }

  data.push_back(w);
  for(int e = 0; e < nAvg; e++){
    data.push_back(quad[e]);
  }
  out.data(data);
}



void scanQuadrupoleMulti(int nAvg, int * l, char * argv0, string name){
  int nThreads = 8;
  vector<thread> thrds;
  OData out(argv0, name);

  for(int i = 0; i < nThreads; i++){
    thrds.push_back(thread(threadCalcQuad, nAvg, nThreads, i, l, ref(out)));
  }

  for(int i = 0; i < nThreads; i++){
    thrds[i].join();
  }

}

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


int main (int arc, char ** argv) {
  //SOTAI
  int l2D1[2] = {10,10};
  int l2D2[2] = {20,20};
  int l2D3[2] = {30,30};
  cout << "Disordered SOTAI" << endl;
  //scanQuadrupoleMulti(40, l2D1, argv[0], "phaseDiagramSOTAI_10x10_m1.1.dat");
  //scanQuadrupoleMulti(40, l2D2, argv[0], "phaseDiagramSOTAI_20x20_m1.1.dat");
  //scanQuadrupoleMulti(1, l2D3, argv[0], "phaseDiagramSOTAI_30x30_m1.1.dat");

  vector<vector<double>> paramList;
  int nPoints = 20;
  for(int i = 0; i <= nPoints; i++){
    vector<double> param; 
    param.push_back(9*(double)i/(double)nPoints);
    paramList.push_back(param);
  }

  /*
     MultiThread r1(quad1, paramList, 8);
     r1.setFile(argv[0], "phaseDiagramSOTAI_10x10_m1.1.dat");
     r1.run();

*/
  MultiThread r2(quad2, paramList, 8);
  r2.setFile(argv[0], "phaseDiagramSOTAI_20x20_m1.1.dat");
  r2.run();
  /*

     MultiThread r3(pol1, paramList, 8);
     r3.setFile(argv[0], "phaseDiagramSOTAIpol_10x10_m1.1.dat");
     r3.run();

     MultiThread r4(pol2, paramList, 8);
     r4.setFile(argv[0], "phaseDiagramSOTAIpol_50x50_m1.1.dat");
     r4.run();
     */

  MultiThread r5(pol3, paramList, 8);
  r5.setFile(argv[0], "phaseDiagramSOTAIpol_100x100_m1.1.dat");
  r5.run();

}
