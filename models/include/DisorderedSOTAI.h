#ifndef DISORDEREDSOTAI_H
#define DISORDEREDSOTAI_H

#include "TBDisorderedH.h"
#include "DisorderFunctions.h"
#include "Wilson.h"

using namespace std;

class DisorderedSOTAI{
  public:
    DisorderedSOTAI(double m, double delta = 0);
    ~DisorderedSOTAI();

    void setOnSite(double);
    void setM(double);

    void setW(double);
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
