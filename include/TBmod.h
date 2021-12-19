#ifndef TBMOD_H
#define TBMOD_H

#include "Hop.h"
#include <armadillo>
#include <vector>

using namespace arma;
using namespace std;

class TBmod{

  public:
    //Tight-binding model with ndim spatial dimensions and norb orbitals per unit cell
    //L is explained in method set_size
    TBmod(int ndim, int norb, int * L = NULL); 
    ~TBmod();

    //Set hopping between norb1 in home cell and orb2 in neighbouring cell defined by vector n
    void set_hop(int norb1, int norb2, int * n, double hop);
    //Set on-site energy for orbital norb
    void set_onsite(int norb, double en);
    //pcb specifies the existance of periodic boundary conditions in each direction (true for PBCs)
    //By default PBCs apply to all spatial directions
    void set_periodic(bool * pbc);
    //Apply twisted boundary conditions with angle (rad) specified by array theta
    //If theta_i = 0 then don't change the boundary conditions in that direction
    void set_twisted(double * theta);
    //Array L specifies the number of unit cells in each direction
    void set_size(int * L);


    //Get Hamiltonian/eigenvalues/eigenvectors in reciprocal space in the directions where PCBs are applied and in real space in the remaining directions
    /*
       cx_mat get_H(double * k);
       vec get_eval(double * k);
       cx_mat get_evec(double * k);
       */
    //Get Hamiltonian in real space with the boundary conditions set by other methods
    /*
       cx_mat get_rH();
       */

  private:

    template<typename T>
      bool check_dim(T * vec);

    vector<Hop> hop;
    int ndim;
    int norb;
    int * L;
    bool * pcb;
    double * theta;

};

#endif
