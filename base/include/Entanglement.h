#ifndef ENTANGLEMENT_H
#define ENTANGLEMENT_H

#include "Hamiltonian.h"

using namespace std;
using namespace arma;

class Entanglement{
  public:
    Entanglement(Hamiltonian * ham, int nOcc, double * k = NULL);
    ~Entanglement();

    void setOcc(int nOcc);

    double bipEntropy(uvec cut);

  private:
    Hamiltonian * ham;
    cx_mat eigStates;
    int nOcc;

    //round numbers close to 0 to 0
    double chop(double);

    void diagonalizeH(double * k = NULL);
};

#endif
