#include <iostream>
#include "DisorderedBBH3D.h"
#include "ParallelMPI.h"

int sampPerJob = 1;
int l[3] = {5,5,5};
double intra = 1.1;

void octu(double * res, double * params){
  DisorderedBBH3D bbh3d(intra, 1);
  bbh3d.setSize(l);
  bbh3d.setW(params[0]);

  for(int i = 0; i < sampPerJob; i++){
    bbh3d.generateDisorder();
    res[i] = bbh3d.getOctupoleManyBody();
  }
}

void quad(double * res, double * params){
  DisorderedBBH3D bbh3d(intra, 1);
  bbh3d.setSize(l);
  bbh3d.setW(params[0]);

  double quadyz;
  double quadxz;
  double quadxy;
  for(int i = 0; i < sampPerJob; i++){
    try{
      bbh3d.generateDisorder();
      quadyz = bbh3d.getBoundQuadrupole(0);
      quadxz = bbh3d.getBoundQuadrupole(1);
      quadxy = bbh3d.getBoundQuadrupole(2);
      res[3*i] = quadyz;
      res[3*i + 1] = quadxz;
      res[3*i + 2] = quadxy;
    }
    catch(const runtime_error & error){
      cout << "Singular matrix found at disorder weight = "  << params[0] << endl;
      i--;
    }
  }
}

int main (int argc, char ** argv) {
  if(argc > 1){
    l[0] = stoi(argv[1]);
    l[1] = l[0];
    l[2] = l[0];
  }

  int sampMult = 40;
  vector<vector<double>> paramList;
  int nPoints = 100;
  for(int i = 0; i <= nPoints; i++){
    vector<double> param; 
    param.push_back(9*(double)i/(double)nPoints);
    paramList.push_back(param);
  }

  if(l[0] == 24){
    paramList.erase(paramList.begin(), paramList.begin() + 95);
  }
  else if(l[0] == 22){
    paramList.erase(paramList.begin(), paramList.begin() + 49);
  }

  ParallelMPI p(&argc, &argv);
  p.setSamples(sampMult);
  p.setParamList(paramList);
  p.setFile(argv[0], "phaseDiagramBBH3Dquad_L" + to_string(l[0]) + "_intra1.1");
  p.setJob(quad, 3*sampPerJob);
  p.run();
}
