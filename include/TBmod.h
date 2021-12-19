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

    //Set hopping between norb1 in home cell and orb2 in neighbouring cell defined by vector n (n should only have non-negative values)
    void set_hop(int norb1, int norb2, int * n, complex<double> hop);
    //Set on-site energy for orbital norb
    void set_onsite(int norb, complex<double> en); //TODO: on site should be separate from hoppings because of h.c.
    //pcb specifies the existance of periodic boundary conditions in each direction (true for PBCs)
    //By default PBCs apply to all spatial directions
    void set_periodic(bool * pbc);
    //Apply twisted boundary conditions with angle (rad) specified by array theta
    //If theta_i = 0 then don't change the boundary conditions in that direction
    void set_twisted(double * theta);
    //Array L specifies the number of unit cells in each direction, L[i] > 0 (if direction i doesn't exist set L[i] = 1)
    void set_size(int * L);


    //Get Hamiltonian/eigenvalues/eigenvectors in reciprocal space in the directions where PCBs are applied and in real space in the remaining directions
    /*
       cx_mat get_H(double * k);
       vec get_eval(double * k);
       cx_mat get_evec(double * k);
       */
    //Get Hamiltonian in real space with the boundary conditions set by other methods
    cx_mat get_rH();

  private:

    //Get Hilbert space index of created orbital in hopping term  
    int get_m1(int * n, Hop & hop);
    //Get Hilbert space index of annihilated orbital in hopping term (PCBs apply)
    int get_m2(int * n, Hop & hop);

    //Calculate accumulated size [Lx,Lx*Ly,Lx*Ly*Lz,...]
    void calc_Laccum();

    //Calculate mesh of unit cells depending on L (full for all directions, otherwise only directions with no PBCs)
    //void calc_n();
    void calc_nfull();

    //Check is system size is big enough for given hopping terms in the directions where PCBs don't apply (hopping can't be larger than system size)
    bool check_size();

    //Hopping and on-site terms
    vector<Hop> hop;
    //System dimension
    int ndim;
    //Number of orbitals per unit cell
    int norb;
    //System size
    int * L;
    //Accumulated system size
    int * Laccum;
    //Periodic boundary conditions
    bool * pcb;
    //Twisted boundary conditions
    double * theta;
    //Mesh of unit cells in all directions
    int ** nfull;
    //Mesh of unit cells only in directions with no PCBs
    int ** n;

};

#endif
