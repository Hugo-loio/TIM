#include <iostream>
#include "DisorderedSOTAI.h"
#include "ParallelMPI.h"

int sampPerJob = 40;
int l[2] = {20,20};
double m = 1.1;

void quad(double * res, double * params){
  DisorderedSOTAI sotai(m);
  sotai.setSize(l);
  sotai.setW(params[0]);
  for(int i = 0; i < sampPerJob; i++){
    sotai.generateDisorder();
    res[i] = sotai.getQuadrupoleManyBody();
  }
}

void pol(double * res, double * params){
  DisorderedSOTAI sotai(m);
  sotai.setSize(l);
  sotai.setW(params[0]);

  for(int i = 0; i < sampPerJob; i++){
    try{
      sotai.generateDisorder();
      double polx = sotai.getBoundPolarization(0);
      double poly = sotai.getBoundPolarization(1);
      res[2*i] = polx;
      res[2*i + 1] = poly;
    }
    catch(const runtime_error & error){
      cout << "Singular matrix found at disorder weight = "  << params[0] << endl;
      i--;
    }
  }
}

int main (int argc, char ** argv) {
  int sampMult = 1;

  vector<vector<double>> paramList;
  int nPoints = 100;
  for(int i = 0; i <= nPoints; i++){
    vector<double> param; 
    param.push_back(9*(double)i/(double)nPoints);
    paramList.push_back(param);
  }

  ParallelMPI p(&argc, &argv);
  p.setSamples(sampMult);
  p.setParamList(paramList);
  //p.setFile(argv[0], "phaseDiagramSOTAI_" + to_string(l[0]) + "x" + to_string(l[1]) + "_m1.1");
  //p.setJob(quad, sampPerJob);
  p.setFile(argv[0], "phaseDiagramSOTAIpol_" + to_string(l[0]) + "x" + to_string(l[1]) + "_m1.1");
  p.setJob(pol, 2*sampPerJob);
  p.run();
}
