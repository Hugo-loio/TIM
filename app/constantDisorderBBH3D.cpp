#include <iostream>
#include "DisorderedBBH3D.h"
#include "ParallelMPI.h"
#include "AuxFunctions.h"

int sampPerJob = 4;
double intra = 1.1;
double weight = 3;

void quad(double * res, double * params){
  DisorderedBBH3D bbh3d(intra, 1);
  int l[3] = {(int)params[0], (int)params[0], (int)params[0]};
  bbh3d.setSize(l);
  bbh3d.setW(weight);

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
      cout << "Singular matrix found" << endl;
      i--;
    }
  }
}

int main (int argc, char ** argv) {
  if(argc > 1){
    weight = stod(argv[1]);
  }
  int version = 0;
  if(argc > 2){
    version = stoi(argv[2]);
  }

  int sampMult = 100;

  int nPoints = 10;
  vector<vector<double>> paramList;
  for(int i = 0; i <= nPoints; i++){
    vector<double> param; 
    param.push_back(20 + i*2);
    paramList.push_back(param);
  }

  for(int i = 0; i < paramList.size(); i++){
    //cout << paramList[i][0] << endl;
  }

  ParallelMPI p(&argc, &argv);
  p.setSamples(sampMult);
  p.setParamList(paramList);
  string fileName = "constantDisorderBBH3Dquad_intra1.1_w" + rmTrailZeros(to_string(weight));
  if(version != 0){
    fileName += "(" + to_string(version) + ")";
  }
  p.setFile(argv[0], fileName + ".dat");
  p.setJob(quad, 3*sampPerJob);
  p.setPrintEachSamp(true);
  p.run();
}
