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
    void setSize(int * l){ham->setSize(l);}
    void setW(double);
    void setLayers(bool *);

    void generateDisorder();

    void getChargeDensity(char * argv0, string fileName, int nx, int ny, int nOrbFilled);
    double getQuadrupoleNestedSupercell(int * l, int * n);
    double getQuadrupoleManyBody();
    double getBoundPolarization(int dir);

    void getSupercellWannierBands(char * argv0, string fileName, int nx, int ny, int dirWilson);

    cx_mat getHam(int * l);

    void test(char * argv0);

  private:
    double m;
    double delta;
    TBModel * model;
    DisorderedHopH2D * ham;
};

#endif
