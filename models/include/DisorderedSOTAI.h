#ifndef DISORDEREDSOTAI_H
#define DISORDEREDSOTAI_H

#include "DisorderedHopH2D.h"

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
    double getQuadrupoleManyBody(int * l);
    double getTopInv(int * l);

    void getSupercellWannierBands(char * argv0, string fileName, int nx, int ny, int dirWilson);

    cx_mat getHam(int * l);

  private:
    double m;
    double delta;
    TBModel * model;
    DisorderedHopH2D * ham;
};

#endif
