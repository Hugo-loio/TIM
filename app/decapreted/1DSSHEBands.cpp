#include <iostream>
#include <armadillo>
#include "SSH_1D.h"
#include "odata.h"

using namespace std;
using namespace arma;

int main(int argc, char ** argv)
{
  int n = 100;
  vector<double> x(n), y1(n), y2(n);

  double t1 = 1;
  double t2 = 2;

  SSH_1D mod(t1,t2);

  for(int i = 0; i < n; i++){
    x[i] = -M_PI+2*M_PI*i/n;
    y1[i] = (mod.eigenvalk(x[i])).at(0);
    y2[i] = (mod.eigenvalk(x[i])).at(1);
  }

  vector<vector<double>> data;
  data.push_back(x);
  data.push_back(y1);
  data.push_back(y2);

  odata od(argv[0], "1DSSHEBands.dat");
  od.data(data);

  return 0;
}

