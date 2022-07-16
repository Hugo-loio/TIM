#include <iostream>
#include <iomanip>
#include <thread>
#include "DisorderedBBH3D.h"
#include "OData.h"
#include "MultiThread.h"
#include "AuxFunctions.h"

void threadQuad(double weight, double gamma, int nSamples, vector<double> & res, vector <double> params){
  DisorderedBBH3D bbh3d(gamma, 1);
  int l[3] = {(int)params[0], (int)params[0], (int)params[0]};
  bbh3d.setSize(l);
  bbh3d.setW(weight);
  res.push_back(params[0]);

  double quadyz;
  double quadxz;
  double quadxy;
  for(int i = 0; i < nSamples; i++){
    try{
      bbh3d.generateDisorder();
      quadyz = bbh3d.getBoundQuadrupole(0);
      quadxz = bbh3d.getBoundQuadrupole(1);
      quadxy = bbh3d.getBoundQuadrupole(2);
      res.push_back(quadyz);
      res.push_back(quadxz);
      res.push_back(quadxy);
    }
    catch(const runtime_error & error){
      cout << "Singular matrix found at disorder weight = "  << params[0] << endl;
    }
  }
}

void quad1(vector<double> & res, vector<double> params){
  threadQuad(3, 1.1, 40, res, params);
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
  for(int i = 0; i <= 10; i++){
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
    cout << paramList[i][0] << endl;
  }

  //run(quad1, paramList, threadNumber, argv[0], "constantDisorderBBH3Dquad_intra1.1_w3", version, part);
}
