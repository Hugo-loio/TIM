#ifndef DISORDEREDBBH3D_H
#define DISORDEREDBBH3D_H

#include "DisorderedHopH3D.h"

using namespace std;

class DisorderedBBH3D{
  public:
    //Intracell hopping t1, intercell hopping t2
    DisorderedBBH3D(double t1 = 1, double t2 = 2, double delta = 0);
    ~DisorderedBBH3D();

    void setOnSite(double);
    void setIntraHop(double);
    void setInterHop(double);
    void setSize(int * l){ham->setSize(l);}
    void setLayers(bool *);
    void setW(double);

    void generateDisorder();

    void getChargeDensity(char * argv0, string fileName, int nOrbFilled);
    double getBoundQuadrupole(int dir);
    double getOctupoleManyBody();

  private:
    TBModel * model;
    DisorderedHopH3D * ham;
};

#endif
