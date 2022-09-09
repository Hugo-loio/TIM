#include "AuxFunctions.h"
#include <iostream>
#include <iomanip>

void run(void (*job) (vector<double> &, vector<double>), vector<vector<double>> & paramList, int threadNumber, char * argv0, string fileName, int nSamples, int version, int part){
  if(part != 0){
    fileName += "_p" + to_string(part);
  }
  if(version != 0){
    fileName +=  "(" + to_string(version) + ")";
  }
  fileName += ".dat";
  MultiThread r(job, paramList, threadNumber);
  r.setFile(argv0, fileName);
  r.setSamples(nSamples);
  r.run();
}

void runSingleThread(void (*job) (vector<double> &, vector<double>), vector<vector<double>> & paramList, char * argv0, string fileName, int nSamples, int version, int part){
  if(part != 0){
    fileName += "_p" + to_string(part);
  }
  if(version != 0){
    fileName +=  "(" + to_string(version) + ")";
  }
  fileName += ".dat";

  MultiThread r(job, paramList, 1);
  r.setFile(argv0, fileName);
  r.setSamples(nSamples);
  r.run();
}

string rmTrailZeros(string str){
  str.erase(str.find_last_not_of('0') + 1, string::npos);
  str.erase(str.find_last_not_of('.') + 1, string::npos);
  return str;
}

double cot(double x){
  return cos(x)/sin(x);
}

void printElapsedTime(chrono::high_resolution_clock::time_point tStart){
  auto tNow = chrono::high_resolution_clock::now();
  auto dH = chrono::duration_cast<chrono::hours>(tNow-tStart);
  auto dM = chrono::duration_cast<chrono::minutes>(tNow-tStart);
  auto dS = chrono::duration_cast<chrono::seconds>(tNow-tStart);
  auto dMS = chrono::duration_cast<chrono::milliseconds>(tNow-tStart);
  int h = dH.count();
  int m = dM.count() % 60;
  int s = dS.count() % 60;
  int ms = dMS.count() % 1000;
  cout << h << ":" << m << ":" << s << ":" << ms << endl;
}
