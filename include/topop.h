#ifndef TOPOP_h
#define TOPOP_h

// Methods for operations on topological matter

#include <armadillo>

namespace topop{

  using namespace arma;

  /*

     Notation for arguments of functions:

   * System with Hamiltonian H(k) in reciprocal space and H in real space
   * Number of discretized k intervals N (loop from k = 0 to k = 2 pi)
   * Number of occupied bands M  
   * Number of discretized twist angle intervals L 

   To get loops with different starting points and directions in reciprocal change argument function H(k) accordingly.
   */

  cx_mat wilson_loop (int N, int M, cx_mat (*H)(double k)); 

  double berry_phase(int N, int M, cx_mat (*H)(double k), double k0); 

  //cx_mat wilson_loop_twisted (int N, int M, int L, cx_mat & H); 

  //double berry_phase_twisted (int N, int M, int L, cx_mat & H);
}

#endif

