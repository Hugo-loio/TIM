#ifndef WILSON_H
#define WILSON_H

#include "Hamiltonian.h"

using namespace std;
using namespace arma;

class Wilson{
  public:
    Wilson(Hamiltonian * ham);
    ~Wilson();

    //n k (momentum) intervals, m occuppied bands
    cx_mat wilsonLoop(int n, int m);
    double berryPhase(int n, int m);
    //n[i]  k (momentum) intervals in direction dir[i], i is order in which the Wilson loops are taken
    //cx_mat nestedWilsonLoop(int * n, int * dir, int m); //TODO: add order

    void setLoopDir(int dir);
    void setLoopStart(double * k0);

    //n theta (twist angle) intervals
    cx_mat wilsonLoopSupercell(int n, int m, double * k = NULL);
    double berryPhaseSupercell(int n, int m, double * k = NULL);

  private:
    double * k0;
    int dir;
    int dim;
    Hamiltonian * ham;
};

#endif
