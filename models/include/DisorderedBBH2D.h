#ifndef DISORDEREDBBH2D_H
#define DISORDEREDBBBH2D_H

#include "TBDisorderedH.h"
#include "DisorderFunctions.h"
#include "Wilson.h"

using namespace std;

class DisorderedBBH2D{
  public:
    //Intracell hopping t1, intercell hopping t2
    DisorderedBBH2D(double t1 = 1, double t2 = 2, double delta = 0);
    ~DisorderedBBH2D();

    void setOnSite(double);
    void setIntraHop(double);
    void setInterHop(double);
    //void setSparse(bool);

    void setProbDisorder(double);
    void generateDisorder();

    void getChargeDensity(char * argv0, string fileName, int nx, int ny, int nOrbFilled);
    double getQuadrupoleNestedSupercell(int * l, int * n);

    void getSupercellWannierBands(char * argv0, string fileName, int nx, int ny, int dirWilson);

  private:
    TBModel * model;
    TBDisorderedH * ham;
    //BoundaryWilsonH * boundH = NULL;
};

#endif
