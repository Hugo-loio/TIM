#include <iostream>
#include <armadillo>
#include "ParallelMPI.h"

using namespace std;
using namespace arma;

void job(double * res, double * params){
  sp_mat mat = sprandu<sp_mat>(1000, 1000, 0.1);
  cx_vec eigVal;
  cx_mat eigVec;

  eigs_gen(eigVal, eigVec, mat, 10, 0.0);
  cout << "MPI test" << endl;
  cout << eigVal << endl;
  res[0] = 0;
}

int main (int argc, char ** argv) {
  int version = 0;
  if(argc > 1){
    version = stoi(argv[1]);
  }

  sp_mat mat = sprandu<sp_mat>(1000, 1000, 0.1);
  cx_vec eigVal;
  cx_mat eigVec;

  //Shift-invert
  eigs_gen(eigVal, eigVec, mat, 10, 0.0);

  cout << eigVal << endl;

  //Test with MPI
  vector<double> param; //Dummy
  param.push_back(0);
  vector<vector<double>> paramList;
  paramList.push_back(param);

  ParallelMPI p(&argc, &argv);
  p.setSamples(2);
  p.setParamList(paramList);
  p.setFile(argv[0], "dummy", version);
  p.setJob(job, 1);
  p.run();
}
