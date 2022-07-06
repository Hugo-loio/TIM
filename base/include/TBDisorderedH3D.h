#ifndef TBDISORDEREDH3D_H
#define TBDISORDEREDH3D_H

#include "TBCleanH.h"

using namespace arma;
using namespace std;

class TBDisorderedH3D : public TBCleanH{
  public:
    TBDisorderedH3D(TBModel model);
    ~TBDisorderedH3D();

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
    complex<double> **** disHop = NULL;
    int nDisOnSite;
    bool * isDisOnSite;
    int * indexDisOnSite;
    complex<double> **** disOnSite = NULL;

    //To be called by generateDisorder of derived class
    void setDisHop(int nDisHop, int * hop);
    void setDisOnSite(int nDisOnSite, int * onSite);

    void createDisArrays();
    void deleteDisArrays();

    template <class mat> void fillDis(mat & res, complex<double> *** w, complex<double> phase, int n, int * incN, int * start, int * end, int dim, int addI = 0, int addJ = 0, int * startI = NULL);
    template <class mat> void fillHop(mat & res, int hopIndex, double * k, int dim, int addI = 0, int addJ = 0, int addN = 0, int * startI = NULL);
    template <class mat> void fillOnSite(mat & res, int onSiteIndex, double *k, int dim, int addN = 0, int * startI = NULL);

    void fillHopWrapper(cx_mat & res, int hopIndex, double *k, int dim, int addI = 0, int addJ = 0, int addN = 0, int * startI = NULL);
    void fillHopWrapper(sp_cx_mat & res, int hopIndex, double *k, int dim, int addI = 0, int addJ = 0, int addN = 0, int * startI = NULL);
    void fillOnSiteWrapper(cx_mat & res, int onSiteIndex, double *k, int dim, int addN = 0, int * startI = NULL);
    void fillOnSiteWrapper(sp_cx_mat & res, int onSiteIndex, double *k, int dim, int addN = 0, int * startI = NULL);
};


#endif
