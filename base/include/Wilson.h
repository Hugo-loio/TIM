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
    vec wilsonPhases(int n, int m);
    cx_mat wilsonEigVec(int n, int m);
    //n[i]  k (momentum) intervals in direction dir[i], dir[i] is order in which the Wilson loops are taken
    cx_mat nestedWilsonLoop(int * n, int * dir, int m); 
    cx_mat nestedWilsonEigVec(int * n, int * dir, int m);
    vec nestedWilsonPhases(int * n, int * dir, int m);
    cx_mat nestedNestedWilsonLoop(int * n, int * dir, int m);

    void setLoopDir(int dir);
    void setLoopStart(double * k0);

    //n theta (twist angle) intervals
    cx_mat wilsonLoopSupercell(int n, int m, double * k = NULL);
    double berryPhaseSupercell(int n, int m, double * k = NULL);
    vec wilsonPhasesSupercell(int n, int m, double * k = NULL);
    cx_mat wilsonEigVecSupercell(int n, int m, double * k = NULL);
    //n[i] theta intervals and m[i] filled states in direction dir[i], dir[i] is order in which the Wilson loops are taken
    cx_mat nestedWilsonLoopSupercell(int * n, int * dir, int * m, double * k = NULL);
    vec nestedWilsonPhasesSupercell(int * n, int * dir, int * m, double * k = NULL);

  private:
    double * k0;
    int dir;
    int dim;
    Hamiltonian * ham;
};

#endif
