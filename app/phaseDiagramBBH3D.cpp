#include <iostream>
#include <iomanip>
#include <thread>
#include "DisorderedBBH3D.h"
#include "OData.h"
#include "MultiThread.h"

void threadOctu(int * l, double gamma, int nSamples, vector<double> & res, vector <double> params){
  DisorderedBBH3D bbh3d(gamma, 1);
  bbh3d.setSize(l);
  bbh3d.setW(params[0]);
  res.push_back(params[0]);

  for(int i = 0; i < nSamples; i++){
    bbh3d.generateDisorder();
    res.push_back(bbh3d.getOctupoleManyBody());
  }
}
void threadQuad(int * l, double gamma, int nSamples, vector<double> & res, vector <double> params){
  DisorderedBBH3D bbh3d(gamma, 1);
  bbh3d.setSize(l);
  bbh3d.setW(params[0]);
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

void run(void (*job) (vector<double> &, vector<double>), vector<vector<double>> & paramList, int threadNumber, char * argv0, string fileName, int version = 0){
  if(version != 0){
    fileName +=  "(" + to_string(version) + ").dat";
  }
  MultiThread r(job, paramList, threadNumber);
  r.setFile(argv0, fileName);
  r.run();
}

void quad1(vector<double> & res, vector<double> params){
  int l[3] = {5,5,5};
  threadQuad(l, 1.1, 40, res, params);
}

void quad2(vector<double> & res, vector<double> params){
  int l[3] = {7,7,7};
  threadQuad(l, 1.1, 40, res, params);
}

void quad3(vector<double> & res, vector<double> params){
  int l[3] = {7,7,7};
  threadQuad(l, 0.9, 40, res, params);
}

void quad4(vector<double> & res, vector<double> params){
  int l[3] = {7,7,7};
  threadQuad(l, 0.5, 40, res, params);
}

void quad5(vector<double> & res, vector<double> params){
  int l[3] = {10,10,10};
  threadQuad(l, 1.1, 40, res, params);
}

void quad6(vector<double> & res, vector<double> params){
  int l[3] = {12,12,12};
  threadQuad(l, 1.1, 40, res, params);
}

void quad7(vector<double> & res, vector<double> params){
  int l[3] = {15,15,15};
  threadQuad(l, 1.1, 40, res, params);
}

int main (int argc, char ** argv) {
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

  //run(quad1, paramList, threadNumber, argv[0], "phaseDiagramBBH3Dquad_L5_intra1.1.dat", version);
  //run(quad2, paramList, threadNumber, argv[0], "phaseDiagramBBH3Dquad_L7_intra1.1.dat", version);
  //run(quad3, paramList, threadNumber, argv[0], "phaseDiagramBBH3Dquad_L7_intra0.9.dat", version);
  //run(quad4, paramList, threadNumber, argv[0], "phaseDiagramBBH3Dquad_L7_intra0.5.dat", version);
  //run(quad5, paramList, threadNumber, argv[0], "phaseDiagramBBH3Dquad_L10_intra1.1.dat", version);
  //run(quad6, paramList, threadNumber, argv[0], "phaseDiagramBBH3Dquad_L12_intra1.1.dat", version);
  //run(quad7, paramList, threadNumber, argv[0], "phaseDiagramBBH3Dquad_L15_intra1.1.dat", version);
}
