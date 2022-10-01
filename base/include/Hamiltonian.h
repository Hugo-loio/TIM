#ifndef HALMILTONIAN_H
#define HALMILTONIAN_H

#include <armadillo>

using namespace arma;

class Hamiltonian{
  public:
    Hamiltonian(int nDim);
    Hamiltonian(const Hamiltonian &);
    ~Hamiltonian();
    virtual cx_mat H(double * k = NULL) = 0;
    virtual sp_cx_mat spH(double * k = NULL) = 0;
    virtual cx_mat blockH(int, int, double * k = NULL) = 0;
    //Does nothing for clean Hamiltonians
    virtual void generateDisorder(){};
    bool getIsSparse(){return isSparse;};
    //nDim twist angle array
    void setTwists(double * theta);
    double * getTwists(){return theta;};
    int getNDim(){return nDim;};
    void setIsReal(bool isReal){this->isReal = isReal;};
    bool getIsReal(){return isReal;};
  protected:
    bool isSparse = true;
    bool isReal = false;
    double * theta;
    int nDim;
};

#endif
