#ifndef HALMILTONIAN_H
#define HALMILTONIAN_H

#include <armadillo>

using namespace arma;

class Hamiltonian{
  public:
    Hamiltonian(int nDim);
    ~Hamiltonian();
    virtual cx_mat H(double * k = NULL) = 0;
    virtual sp_cx_mat spH(double * k = NULL) = 0;
    bool getIsSparse(){return isSparse;};
    //nDim twist angle array
    void setTwists(double * theta);
  protected:
    bool isSparse = true;
    double * theta;
    int nDim;
};

#endif
