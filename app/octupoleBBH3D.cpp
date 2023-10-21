#include <iostream>
#include "DisorderedBBH3D.h"
#include "ParallelMPI.h"
#include "AuxFunctions.h"

int sampPerJob = 1;
int l[3] = {5,5,5};
double intra = 1.1;
double maxW = 9;

void octu(double * res, double * params){
  DisorderedBBH3D bbh3d(intra, 1);
  bbh3d.setSize(l);
  bbh3d.setW(params[0]);

  for(int i = 0; i < sampPerJob; i++){
    bbh3d.generateDisorder();
    res[i] = bbh3d.getOctupoleManyBody();
  }
}


int main (int argc, char ** argv) {
  int version = 0;
  int sampMult = 10;

  if(argc > 3){
    l[0] = stoi(argv[1]);
    l[1] = l[0];
    l[2] = l[0];
    sampMult = stoi(argv[2]);
    version = stoi(argv[3]);
  }

  vector<vector<double>> paramList;
  //vector<double> param; 
  //param.push_back(maxW*(double)i/(double)nPoints);
  paramList.push_back(vector<double>(1, 1));
  paramList.push_back(vector<double>(1, 2));
  paramList.push_back(vector<double>(1, 2.7));
  paramList.push_back(vector<double>(1, 3));
  paramList.push_back(vector<double>(1, 3.2));

  ParallelMPI p(&argc, &argv);
  p.setSamples(sampMult);
  p.setParamList(paramList);
  p.setFile(argv[0], "octupoleBBH3D_L" + to_string(l[0]) + "_intra" + rmTrailZeros(to_string(intra)), version);
  p.setJob(octu, sampPerJob);
  p.run();
}
