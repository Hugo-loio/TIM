#ifndef TBCLEANH_H
#define TBCLEANH_H

#include "TBModel.h"
#include "Hamiltonian.h"

using namespace arma;
using namespace std;

class TBCleanH : public Hamiltonian{
  public:
    TBCleanH(TBModel model);
    TBCleanH(const TBCleanH &);
    ~TBCleanH();

    //Set boundary conditions for each direction: 0 for open, 1 for periodic (goes to reciprocal space), 2 for twisted
    void setBC(int * bC);
    //nDim array that specifies the number of unit cells in each direction (l[i] > 0)    
    void setSize(int * l);
    const int * getSize(){return l;}

    //Get Hamiltonian in reciprocal space in the directions where PCBs are applied and in real space in the remaining directions
    cx_mat H(double * k = NULL);
    sp_cx_mat spH(double * k = NULL);
    //Using open boundary conditions for the directions not included in the blocks
    cx_mat blockH(int line, int col, double * k = NULL);

    //Whether the user prefers Sparse matrices or not
    void setSparse(bool);

    //Order of the spatial directions in the Hamiltonian matrix
    void setOrder(int * o);

    //Orbital layer coordinate within unit cell
    void setOrbLayer(int ** l);

    //Dimensionality  of the layers associated with each block
    void setBlockDim(int bDim);

    int getIndex(int orb, int * n){
      if(!isUpdated){
	calcAux();
	isUpdated = true;
      }

      int res = flatten(orb, n);
      return res;
    };

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
    //Dimensions  of the layers in each block (nRDim > bDim) 
    int bDim;

    //Auxiliary pre-calculations (boost performance)
    int * nHopBulk;
    int ** nHopBound;
    int * inc;
    int ** incNHopBulk;
    int *** incNHopBound;
    // For lattice with sites
    int ** startHopBulk;
    int ** startHopBound;
    int ** endHopBulk;
    int ** endHopBound;
    // For unit cells
    int ** startHopUCBulk;
    int ** startHopUCBound;
    int ** endHopUCBulk;
    int ** endHopUCBound;
    int ** incNOnSite;
    int ** startOnSite;
    int ** endOnSite;
    //Index of orbital inside cell subdivision 
    int * mOrb;
    //Number of layers per unit cell in each direction
    int * nu;
    //Whether the auxiliary pre-calculations have been preformed for the current options
    bool isUpdated = false;
    //Number of layers multiplied (definition in notes)
    int * nuAccum;
    //Accumulated system size 
    int * lAccum; 

    //Calculate auxiliary quantities
    void calcAux();
    //Get flattened orbital index
    int flatten(int alpha, int * n);
    // for already sorted n
    int flatten2(int alpha, int * n);

    template <class mat> void fill(mat & res, complex<double> w, int n, int * incN, int * start, int * end, int dim, int addI = 0, int addJ = 0, int * startI = NULL);
    template <class mat> void fillHop(mat & res, int hopIndex, double * k, int dim, int addI = 0, int addJ = 0, int addN = 0, int * startI = NULL);
    template <class mat> void fillOnSite(mat & res, int onSiteIndex, double *k, int dim, int addN = 0, int * startI = NULL);

    //Wrappers to override in derived class
    virtual void fillHopWrapper(cx_mat & res, int hopIndex, double *k, int dim, int addI = 0, int addJ = 0, int addN = 0, int * startI = NULL);
    virtual void fillHopWrapper(sp_cx_mat & res, int hopIndex, double *k, int dim, int addI = 0, int addJ = 0, int addN = 0, int * startI = NULL);
    virtual void fillOnSiteWrapper(cx_mat & res, int onSiteIndex, double *k, int dim, int addN = 0, int * startI = NULL);
    virtual void fillOnSiteWrapper(sp_cx_mat & res, int onSiteIndex, double *k, int dim, int addN = 0, int * startI = NULL);
};

#endif
