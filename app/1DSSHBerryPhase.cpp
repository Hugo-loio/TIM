#include <iostream>
#include "matplotlibcpp.h"
#include <armadillo>
#include "SSH_1D.h"

using namespace std;
using namespace arma;
namespace plt = matplotlibcpp;

int main()
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

  plt::plot(x,y);
  plt::xlabel("J'");
  plt::ylabel("$ \\phi$");
  plt::title("Berry phase, J = 1");
  plt::show();

  return 0;
}

