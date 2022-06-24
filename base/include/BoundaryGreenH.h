#ifndef BOUNDARYGREENH_H
#define BOUNDARYGREENH_H

#include "Hamiltonian.h"

using namespace arma;
using namespace std;

class BoundaryGreenH : public Hamiltonian{
  public:
    BoundaryGreenH(Hamiltonian * ham, int blockSize, int nLayers);
    ~BoundaryGreenH();

    cx_mat boundaryGreenFunc(double * k = NULL);
    cx_mat H(double * k = NULL);
    //Dummy sparse option
    sp_cx_mat spH(double * k = NULL){return sp_cx_mat(H(k));};
    cx_mat blockH(int line, int col, double * k = NULL){return cx_mat();};
    void setLowerBound(bool lowerBound, int blockSize);

  private:
    bool lowerBound = false;
    int lowerBlockSize;
    int blockSize;
    int nLayers;
    Hamiltonian * ham;
};

#endif
