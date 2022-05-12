#ifndef MULTIPOLEOP_H
#define MULTIPOLEOP_H

#include "Hamiltonian.h"

using namespace std;
using namespace arma;

class MultipoleOp{
  public:
    //System size l[i] in direction i, nOrb atomic orbitals per unit cell, dim dimensions in real space
    MultipoleOp(Hamiltonian * ham, int * l, int dim, int nOrb = 1);
    ~MultipoleOp();

    void setOcc(int nOcc){this->nOcc = nOcc;};
    //Directions a,b,c
    double polarization(int a, double * k = NULL);
    double quadrupole(int a, int b, double * k = NULL);
    double octupole(int a, int b, int c, double * k = NULL);

  private:
    Hamiltonian * ham;
    int * lAccum;
    int * l;
    int dim;
    int nOcc;
    int nOrb;


    //Get integer index corresponding to unit cell number vector
    int getN(int * n);
    //Recursive function that takes point and changes it to next point in mesh
    void nextPoint(int depth, int * point, bool up);
};

#endif
