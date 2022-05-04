#ifndef DISORDEREDBBH3D_H
#define DISORDEREDBBH3D_H

#include "TBDisorderedH.h"
#include "DisorderFunctions.h"
#include "Wilson.h"

using namespace std;

class DisorderedBBH3D{
  public:
    //Intracell hopping t1, intercell hopping t2
    DisorderedBBH3D(double t1 = 1, double t2 = 2, double delta = 0);
    ~DisorderedBBH3D();

    void setOnSite(double);
    void setIntraHop(double);
    void setInterHop(double);
    //void setSparse(bool);

    void setProbDisorder(double);
    void generateDisorder();

    void getChargeDensity(char * argv0, string fileName, int nx, int ny, int nOrbFilled);

    double getOctupoleNestedSupercell(int * l, int * n);
    void getSupercellNestedWannierBands(char * argv0, string fileName, int * l, int * n);

  private:
    TBModel * model;
    TBDisorderedH * ham;
    //BoundaryWilsonH * boundH = NULL;
};

#endif
