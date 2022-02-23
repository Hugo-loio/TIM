#ifndef TBCLEANH_H
#define TBCLEANH_H

#include "TBModel.h"
#include "Hamiltonian.h"

using namespace arma;
using namespace std;

class TBCleanH : public Hamiltonian{
  public:
    TBCleanH(TBModel & model);
    //TBCleanH(const TBCleanH &);
    ~TBCleanH();

    //Set boundary conditions for each direction: 0 for open, 1 for periodic (goes to reciprocal space), 2 for twisted
    void setBC(int * bC);
    //nDim twist angle array used in directions with twisted boundary conditions
    void setTwists(double * theta);
    //nDim array that specifies the number of unit cells in each direction (l[i] > 0)    
    void setSize(int * l);

    //Get Hamiltonian in reciprocal space in the directions where PCBs are applied and in real space in the remaining directions
    cx_mat H(double * k);
    //sp_cx_mat spH(double * k);

    //Whether the user prefers Sparse matrices or not
    void setSparse(bool);

  private:
    //Information on TB model
    TBModel model;
    int nDim;
    int nOrb;
    int nHop;
    int nOnSite;

    //Number of directions with no PBCs
    int nRDim;
    //Index of directions with no PBCs
    int * rIndex;
    //System size
    int * l;
    //Size of boundary (depends on the max neighbour order of the hopping terms)
    int ** lBound;
    //Accumulated system size [Lx, LxLy, ...] (only in directions with no PCBs)
    int * lAccum; 
    //Boundary volume 
    int vBound;
    //Bulk volume
    int vBulk;
    //Boundary conditions
    int * bC;
    //Twists for twisted boundary conditions
    double * theta;
    //Meshes of unit cells with nDim spatial indexes and one collapsed 1D index 
    int ** nBulk;
    int ** nBound;

    //Get integer index corresponding to unit cell number vector
    int getN(int * n);
    //Generate meshes of unit cells
    void calcN();
    //Recursive function that takes point and changes it to next point in mesh
    void nextPoint(int depth, int * point, bool up);
    //Check is system size is big enough for given hopping terms (increase size if it isn't)
    void checkSize();
    //Calculate volumes
    void calcVol();
};

#endif
