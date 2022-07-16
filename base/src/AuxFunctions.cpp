#include "AuxFunctions.h"

void run(void (*job) (vector<double> &, vector<double>), vector<vector<double>> & paramList, int threadNumber, char * argv0, string fileName, int version, int part){
  if(part != 0){
    fileName += "_p" + to_string(part);
  }
  if(version != 0){
    fileName +=  "(" + to_string(version) + ").dat";
  }
  MultiThread r(job, paramList, threadNumber);
  r.setFile(argv0, fileName);
  r.run();
}

