#ifndef WILSON_H
#define WILSON_H

#include "Hamiltonian.h"

using namespace std;
using namespace arma;

class Wilson{
  public:
    Wilson(Hamiltonian * ham, int dim);
    ~Wilson();

    //n k (momentum) intervals, m occuppied bands
    cx_mat wilsonLoop(int n, int m);
    double berryPhase(int n, int m);

    void setLoopDir(int dir);
    void setLoopStart(double * k0);

    //n theta (twist angle) intervals
    //cx_mat wilsonLoopSupercell(int n, int m);
    //double berryPhaseSupercell(int n, int m);

  private:
    double * k0;
    int dir;
    int dim;
    Hamiltonian * ham;

    cx_mat evOcc(int m, cx_mat h);
};

#endif
