#ifndef DISORDEREDSOTAI_H
#define DISORDEREDSOTAI_H

#include "DisorderedHopH2D.h"
#include "LocalizationProps.h"
#include "DOS.h"

using namespace std;

class DisorderedSOTAI{
  public:
    DisorderedSOTAI(double m, double delta = 0);
    ~DisorderedSOTAI();

    void setOnSite(double);
    void setM(double);
    void setSize(int * l);
    void setW(double);
    void setLayers(bool *);

    void generateDisorder();

    void getChargeDensity(char * argv0, string fileName, int nx, int ny, int nOrbFilled);
    double getQuadrupoleNestedSupercell(int * l, int * n);
    double getQuadrupoleManyBody();
    double getBoundPolarization(int dir);

    double getIPR(int nStates, double en = 0);
    vector <double> getTMM(int qrIt, double en, int m, int dir);
    double getLSR(int nStates, double en = 0);

    double getDOS(double en, int nMoments, int nRandVecs, double eMax = 0);
    double getLDOS(int * n, double en, int nMoments, double eMax = 0);
    double getLDOS(int * n, double en, double range);

    void getSupercellWannierBands(char * argv0, string fileName, int nx, int ny, int dirWilson);

    double getEnGap(double en);

    double getMaxE();

    cx_mat getHam();
    void printHam(char * argv0, string fileName);

    void test(char * argv0);

  private:
    TBModel * model;
    DisorderedHopH2D * ham;
    DOS * dos = NULL;
    LocalizationProps * loc;

    bool updateDOS = true;
};

#endif
