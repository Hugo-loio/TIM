#ifndef TBMODOP_H
#define TBMODOP_H

#include "TBmod.h"

using namespace arma;
using namespace std;

/*
   Notation for arguments of functions:

 * System with Hamiltonian H(k) in reciprocal space and H in real space
 * Number of discretized k intervals N (loop from k = 0 to k = 2 pi)
 * Number of occupied bands M  
 * Number of discretized twist angle intervals L 

*/

class TBmodOp: public TBmod{
  public:
    TBmodOp(int ndim, int norb, int * L = NULL);
    //TODO: Copy const
    ~TBmodOp();

    //TODO: Generalized kloop (with some interpolation)
    cx_mat wilson_loop_k(int N, int M);
    double berry_phase_k(int N, int M);

    //cx_mat wilson_loop_twisted(int N, int M, int L);
    //double berry_phase_twisted(int N, int M, int L);

    //Use sparse or regular Hamiltonian matrices
    void set_sparse(bool);
    //Direction in which quantities are calculated
    void set_dir(int dir);
    //k at start of loop
    void set_k0(double * k0);

  private:

    //Se if returning by reference makes diference in speed
    sp_cx_mat sp_hloop(double k);
    cx_mat hloop(double k);

    //cx_mat htwist(double theta);
    //sp_cx_mat sp_htwist(double theta);

    cx_mat ev_occ(int M, cx_mat H);

    bool sparse;
    int dir;
    double * k0;
};

#endif
