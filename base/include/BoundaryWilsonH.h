#ifndef BOUNDARYWILSONH_H
#define BOUNDARYWILSONH_H

#include "Hamiltonian.h"
#include "Wilson.h"

class BoundaryWilsonH : public Hamiltonian{
  public:
    //n loop steps, m occupied bands
    BoundaryWilsonH(Hamiltonian * h, int dir, int n, int m);
    ~BoundaryWilsonH();

    cx_mat H(double * k = NULL);
    sp_cx_mat spH(double * k = NULL);

  private:
    //Parent hamiltonian
    Hamiltonian * ham;
    Wilson wilson;
    int dir;
    int nOcc;
    int nK;
};

#endif
