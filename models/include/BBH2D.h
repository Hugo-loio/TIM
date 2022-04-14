#ifndef BBH2D_H
#define BBH2D_H

#include "TBCleanH.h"
#include <armadillo>
#include <string>
#include "Wilson.h"

using namespace arma;
using namespace std;

class BBH2D{
  public:
    //Intracell hopping t1, intercell hopping t2
    BBH2D(double t1 = 1, double t2 = 2, double delta = 0);
    ~BBH2D();

    void setOnSite(double);
    void setIntraHop(double);
    void setInterHop(double);
    //void setSparse(bool);

    double berryPhase(int n, int dir = 0, double * k0 = NULL);
    void getBands(char * argv0, string fileName, int nx, int ny);
    void getWannierBands(char * argv0, string fileName, int dir);
    void getChargeDensity(char * argv0, string fileName, int nx, int ny, int nOrbFilled);
    //cx_mat getH(double * k = NULL);
    double getQuadrupoleNested(int nx, int ny, double * k0 = NULL);

    void getSupercellWannierBands(char * argv0, string fileName, int nx, int ny, int dirWilson);

  private:
    TBModel * model;
    TBCleanH * ham;
    //BoundaryWilsonH * boundH = NULL;
};

#endif
