#ifndef BOUNDARYGREENH_H
#define BOUNDARYGREENH_H

#include "Hamiltonian.h"

using namespace arma;
using namespace std;

class BoundaryGreenH : public Hamiltonian{
  public:
    BoundaryGreenH(Hamiltonian * ham, int nOrb, int * l);
    ~BoundaryGreenH();

    cx_mat boundaryGreenFunc(double * k);
    cx_mat H(double * k);
    //Dummy sparse option
    sp_cx_mat spH(double * k){return sp_cx_mat();};

  private:
    int blockSize;
    int nLayers;
    Hamiltonian * ham;
};

#endif
