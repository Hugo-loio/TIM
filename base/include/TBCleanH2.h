#ifndef TBCLEANH2_H
#define TBCLEANH2_H

#include "TBModel.h"
#include "Hamiltonian.h"

using namespace arma;
using namespace std;

class TBCleanH2 : public Hamiltonian{
  public:
    TBCleanH2(TBModel model);
    //TBCleanH(const TBCleanH &);
    ~TBCleanH2();

    //Set boundary conditions for each direction: 0 for open, 1 for periodic (goes to reciprocal space), 2 for twisted
    void setBC(int * bC);
    //nDim array that specifies the number of unit cells in each direction (l[i] > 0)    
    void setSize(int * l);

    //Get Hamiltonian in reciprocal space in the directions where PCBs are applied and in real space in the remaining directions
    cx_mat H(double * k = NULL);
    sp_cx_mat spH(double * k = NULL);

    //Whether the user prefers Sparse matrices or not
    void setSparse(bool);

    //Order of the spatial directions in the Hamiltonian matrix
    void setOrder(int * o);

    //Orbital layer coordinate within unit cell
    void setOrbLayer(int ** l);

  protected:
    //Information on TB model
    TBModel model;
    int nOrb;
    int nHop;
    int nOnSite;

    //Number of directions with no PBCs
    int nRDim;
    //Order of directions
    int * rOrder;
    //Layer coordinates of each orbital
    int ** lOrb;
    //Index of ordered directions with no PBCs
    int * rIndex;
    //System size
    int * l;
    //Boundary conditions
    int * bC;

    //Auxilary pre-calculations (boost performance)
    int * nHopBulk;
    int ** nHopBound;
    int ** incHop;
    int ** incNHop;
    int ** startHopBulk;
    int ** startHopBound;
    int ** endHopBulk;
    int ** endHopBound;
    int ** incOnSite;
    int ** incNOnSite;
    int ** startOnSite;
    int ** endOnSite;
    //Index of orbital inside cell subdivision 
    int * mOrb;
    //Number of layers per unit cell in each direction
    int * nu;
    //Whether the auxiliary pre-calculations have been preformed for the current options
    bool isUpdated;
    //Number of layers multiplied
    int * nuAccum;
    //Accumulated system size 
    int * lAccum; 

    //Calculate auxiliary quantities
    void calcAux();
    //Get flattened orbital index
    int flatten(int alpha, int * n);
};

#endif
