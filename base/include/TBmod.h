#ifndef TBMOD_H
#define TBMOD_H

#include "Hop.h"
#include "Onsite.h"
#include <armadillo>
#include <vector>

using namespace arma;
using namespace std;

class TBmod{

  public:
    //Tight-binding model with ndim spatial dimensions and norb orbitals per unit cell
    //L is explained in method set_size
    TBmod(int ndim, int norb, int * L = NULL); 
    //TODO: copy const
    ~TBmod();

    //Set hopping between norb1 in home cell and orb2 in neighbouring cell defined by vector n 
    void set_hop(int norb1, int norb2, int * n, complex<double> hop);
    //Set on-site energy for orbital norb
    void set_onsite(int norb, complex<double> en); 
    //Set boundary conditions for each direction: 0 for open, 1 for periodic (default BCs), 2 for twisted
    void set_bc(int * bc);
    //ndim twist angle array used in directions with twisted boundary conditions
    void set_twists(double * theta);
    //Array L specifies the number of unit cells in each direction, L[i] > 0 (if direction i doesn't exist set L[i] = 1)
    void set_size(int * L);
    //Sparse or normal matrix to be used for Hamiltonian
    //void set_sparse(bool val);


    //Get Hamiltonian/eigenvalues/eigenvectors in reciprocal space in the directions where PCBs are applied and in real space in the remaining directions
    cx_mat get_H(double * k);
    sp_cx_mat get_spH(double * k);
    /*
       vec get_eval(double * k);
       cx_mat get_evec(double * k);
       */
    //Get Hamiltonian in real space with the boundary conditions set by other methods
    cx_mat get_rH();
    sp_cx_mat get_sprH();


  private:

    //Get integer index corresponding to unit cell position vector
    int get_n(int * n);
    int get_nfull(int * n);
    //Generate meshes of unit cells
    void calc_n();
    void calc_nfull();
    //Recursive functions that take point and change it to next point in mesh
    void next_point_n(int depth, int * point, bool up);
    void next_point_nfull(int depth, int * point, bool up);
    //Check is system size is big enough for given hopping terms (increase size if it isn't)
    void check_size();
    //Calculate volumes
    void calc_vol();

    //Hopping and on-site terms
    vector<Hop> hop;
    vector<Onsite> os; 
    //Sparse matrix for hamiltonian

    //System dimension
    int ndim;
    //Number of directions with no PBCs
    int nrdim;
    //Index of directions with no PBCs
    int * rindex;
    //Number of orbitals per unit cell
    int norb;
    //System size
    int * L;
    //Size of boundary (depends on the max neighbour order of the hopping terms)
    int ** Lbound;
    //Accumulated system size [Lx, LxLy, ...]
    int * Laccum; 
    //Accumulated system size only in directions with no PBCs
    int * Lraccum;
    //Boundary volume 
    int Vbound;
    int Vboundfull;
    //Bulk volume
    int Vbulk;
    int Vbulkfull;
    //Boundary conditions
    int * bc;
    //Twists for twisted boundary conditions
    double * theta;
    //Meshes of unit cells with ndim spatial indexes and one collapsed 1D index 
    //Mesh of unit cells in all directions (bulk and boundary)
    int ** nfull_bulk;
    int ** nfull_bound;
    //Mesh of unit cells only in directions with no PCBs (bulk and boundary)
    int ** n_bulk;
    int ** n_bound;
};

#endif
