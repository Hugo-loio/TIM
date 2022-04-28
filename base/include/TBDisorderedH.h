#ifndef TBDISORDEREDH_H
#define TBDISORDEREDH_H

#include "TBCleanH.h"

using namespace arma;
using namespace std;

class TBDisorderedH : public TBCleanH{
  public:
    TBDisorderedH(TBModel model);
    ~TBDisorderedH();

    //Set boundary conditions for each direction: 0 for open, 1 for periodic (goes to reciprocal space), 2 for twisted
    void setBC(int * bC);
    //nDim array that specifies the number of unit cells in each direction (l[i] > 0)    
    void setSize(int * l);

    //Get Hamiltonian in reciprocal space in the directions where PCBs are applied and in real space in the remaining directions
    cx_mat H(double * k);
    sp_cx_mat spH(double * k);

    //func receives Hop/OnSite object and outputs a disordered value for that hop/on-site term
    //If func == NULL, the hop/on-site terms will correspond to the clean system
    void setHopDisorderFunction(complex<double> (*func)(Hop &, double) = NULL){hopFunc = func;};
    void setOnSiteDisorder(complex<double> (*func)(OnSite &, double) = NULL){onSiteFunc = func;};

    //Set weights for disorder functions
    void setHopDisorderWeight(double w){wHop = w;}
    void setOnSiteDisorderWeight(double w){wOnSite = w;}

    void generateDisorder();
  private:

    void resize();

    double wHop;
    double wOnSite;

    complex<double> (* hopFunc)(Hop &, double) = NULL;
    complex<double> (* onSiteFunc)(OnSite &, double) = NULL;

    complex<double> ** hopBulk;
    complex<double> ** hopBound;
    complex<double> ** onSiteBulk;
    complex<double> ** onSiteBound;

    void delete_arrays();
};

#endif
