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
  vector<double> x(n), y1(n), y2(n);

  double t1 = 1;
  double t2 = 2;

  SSH_1D mod(t1,t2);

  for(int i = 0; i < n; i++){
    x[i] = -M_PI+2*M_PI*i/n;
    y1[i] = (mod.eigenvalk(x[i])).at(0);
    y2[i] = (mod.eigenvalk(x[i])).at(1);
  }

  plt::plot(x,y1);
  plt::plot(x,y2);
  plt::xlabel("k");
  plt::ylabel("E");
  plt::title("1D SSH energy bands, J = 1, J' = 2");
  plt::show();

  return 0;
}

