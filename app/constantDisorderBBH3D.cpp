#include <iostream>
#include "DisorderedBBH3D.h"
#include "ParallelMPI.h"

int sampPerJob = 10;
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
  int part = 0;
  if(argc > 1){
    part = stoi(argv[1]);
  }

  int sampMult = 10;

  int nPoints = 10;
  vector<vector<double>> paramList;
  for(int i = 0; i <= nPoints; i++){
    vector<double> param; 
    param.push_back(10 + i*2);
    paramList.push_back(param);
  }

  int nParts = 4;
  if(part != 0){
    switch(part){
      case 1:
	paramList.erase(paramList.begin() + 5, paramList.end());
	break;
      case 2:
	paramList.erase(paramList.begin(), paramList.begin() + 5);
	paramList.erase(paramList.begin() + 3, paramList.end());
	break;
      case 3:
	paramList.erase(paramList.begin(), paramList.begin() + 8);
	paramList.erase(paramList.begin() + 2, paramList.end());
	break;
      case 4:
	paramList.erase(paramList.begin(), paramList.begin() + 10);
	break;
    }
  }

  for(int i = 0; i < paramList.size(); i++){
    //cout << paramList[i][0] << endl;
  }

  ParallelMPI p(&argc, &argv);
  p.setSamples(sampMult);
  p.setParamList(paramList);
  string fileName = "constantDisorderBBH3Dquad_intra1.1_w3";
  if(part != 0){
    fileName += "p" + to_string(part);
  }
  p.setFile(argv[0], fileName + ".dat");
  p.setJob(quad, 3*sampPerJob);
  p.run();
}
