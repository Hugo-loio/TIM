#include <iostream>
#include <armadillo>
#include "SSH_1D.h"
#include "odata.h"

using namespace std;
using namespace arma;

int main(int argc, char * argv[])
{
  int n = 100;
  int N = 10;
  vector<double> x(n), y(n);

  double t1 = 1;
  double t2 = 2;

  SSH_1D mod(t1,t2);

  for(int i = 0; i < n; i++){
    x[i] = 2*(double)i/n;
    mod.set_interhop(x[i]);
    y[i] = mod.berry_kspace(N); 
  }

  vector<vector<double>> data;
  data.push_back(x);
  data.push_back(y);
  odata od(argv[0], "1DSSHBerryPhase.dat");
  od.data(data);

  cout << argv[0] << endl;

  return 0;
}

