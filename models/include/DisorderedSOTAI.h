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

    double getIPR(int nStates);
    double getTMM(int nIt, double en, int l);

    void getSupercellWannierBands(char * argv0, string fileName, int nx, int ny, int dirWilson);

    cx_mat getHam();
    void printHam(char * argv0, string fileName);

    void test(char * argv0);

  private:
    TBModel * model;
    DisorderedHopH2D * ham;
};

#endif
