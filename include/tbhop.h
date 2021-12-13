#ifndef TBHOP_h
#define TBHOP_h

// Methods for operations on tight binding hamiltonians

#include <armadillo>

namespace tbhop{

  using namespace arma;

  /*
     Notation for arguments of functions:

   * System with Hamiltonian H(k) in reciprocal space and H in real space
   * Number of discretized k intervals N (loop from k = 0 to k = 2 pi)
   * Number of occupied bands M  
   * Number of discretized twist angle intervals L 

*/

  // Returns matrix in which the columns are the eigenvectors of the occuppied states
  cx_mat ev_occ(int M, cx_mat H);


  // Abstract classes which models must inherit from in order to use some of the methods 
  class Kloop{

    protected:
      //Hamiltonian defined along a loop in the BZ  parameterized by k (0 to 2 pi)
      virtual cx_mat hloop(double k) = 0;

      // Wilson loop over path defined by h
      cx_mat wilson_loop (int N, int M); 

      // Berry phase over path defined by h
      double berry_phase(int N, int M); 
  };


  class twisted{
    public:
      //Hamiltonian defined in real space with twisted boundary conditions of angle theta
      virtual cx_mat htwist(double theta) = 0;

      //cx_mat wilson_loop_twisted (int N, int M, int L); 

      //double berry_phase_twisted (int N, int M, int L);
  };


  // TODO: Function that checks validity of arguments?

}

#endif

