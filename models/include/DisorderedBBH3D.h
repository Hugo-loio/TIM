#ifndef DISORDEREDBBH3D_H
#define DISORDEREDBBH3D_H

#include "DisorderedHopH3D.h"
#include "LocalizationProps.h"
#include "DOS.h"

using namespace std;

class DisorderedBBH3D{
  public:
    //Intracell hopping t1, intercell hopping t2
    DisorderedBBH3D(double t1 = 0.5, double t2 = 1, double delta = 0);
    ~DisorderedBBH3D();

    void setOnSite(double);
    void setIntraHop(double);
    void setInterHop(double);
    void setSize(int * l);
    void setLayers(bool *);
    void setW(double);

    void generateDisorder();

    void getChargeDensity(char * argv0, string fileName, int nOrbFilled);
    double getBoundQuadrupole(int dir);
    double getOctupoleManyBody();
    double getIPR(int nStates, double en = 0);
    double getTMM(int qrIt, double en, int m);
    double getLSR(int nStates, double en = 0);
    double getDOS(double en, int nMoments, int nRandVecs, double eMax = 0);
    double getEnGap(double en);
    double getMaxE();

  private:
    TBModel * model;
    DisorderedHopH3D * ham;
    DOS * dos = NULL;
    LocalizationProps * loc;

    bool updateDOS = true;
};

#endif
