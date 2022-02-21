#ifndef HALMILTONIAN_H
#define HALMILTONIAN_H

#include <armadillo>

using namespace arma;

class Hamiltonian{
  public:
    virtual cx_mat H(double * k = NULL) = 0;
    virtual sp_cx_mat spH(double * k = NULL) = 0;
  private:
    bool isSparse = true;
};

#endif
