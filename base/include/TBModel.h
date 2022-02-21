#ifndef TBMODEL_H
#define TBMODEL_H

#include "Hop.h"
#include "OnSite.h"
#include <armadillo>
#include <vector>

using namespace arma;
using namespace std;

class TBModel{
  public:
    //Tight-binding model with nDim spatial dimensions and nOrb orbitals per unit cell
    TBModel(int nDim, int nOrb); 
    TBModel(const TBModel &);
    ~TBModel();

    void setHop(int nOrb1, int nOrb2, int * n, complex<double> hop);
    void setHop(Hop hop);
    void setOnSite(int nOrb, complex<double> en); 
    void setOnSite(OnSite onSite); 
    Hop & getHop(int nHop);
    OnSite & getOnSite(int nOnSite);

  protected:

    //Hopping and on-site terms
    vector<Hop> hop;
    vector<OnSite> onSite; 

    //System dimension
    int nDim;
    //Number of orbitals per unit cell
    int nOrb;
};

#endif
