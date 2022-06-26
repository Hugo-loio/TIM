#ifndef TBDISORDEREDH2D
#define TBDISORDEREDH2D

#include "TBCleanH.h"

using namespace arma;
using namespace std;

class TBDisorderedH2D : public TBCleanH{
  public:
    TBDisorderedH2D(TBModel model);
    ~TBDisorderedH2D();

    void setSize(int * l);

    //Fill disHop and disOnSite
    virtual void generateDisorder() = 0;

    cx_mat H(double * k = NULL);
    sp_cx_mat spH(double * k = NULL);
    cx_mat blockH(int line, int col, double * k = NULL);

  protected:

    int nDisHop;
    int * indexDisHop;
    bool * isDisHop;
    complex<double> *** disHop = NULL;
    int nDisOnSite;
    bool * isDisOnSite;
    int * indexDisOnSite;
    complex<double> *** disOnSite = NULL;

    //To be called by generateDisorder of derived class
    void setDisHop(int nDisHop, int * hop);
    void setDisOnSite(int nDisOnSite, int * onSite);

    void createDisArrays();
    void deleteDisArrays();
};

#endif
