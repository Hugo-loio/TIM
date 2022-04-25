#ifndef BBH3D_H
#define BBH3D_H

#include "TBCleanH.h"
#include <armadillo>
#include <string>
#include "Wilson.h"

using namespace arma;
using namespace std;

class BBH3D{
  public:
    //Intracell hopping t1, intercell hopping t2
    BBH3D(double t1 = 1, double t2 = 2, double delta = 0);
    ~BBH3D();

    void setIntraHop(double);
    void setInterHop(double);
    void setOnSite(double);
    //void setSparse(bool);

    //double berryPhase(int n, int dir = 0, double * k0 = NULL);
    void getBands(char * argv0, string fileName);
    void getChargeDensity(char * argv0, string fileName, int * l, int nOrbFilled);
    void getWannierBands(char * argv0, string fileName, int dir);
    void getNestedWannierBands(char * argv0, string fileName);
    cx_mat getH(double * k = NULL);
    double getOctupoleNested(int nx, int ny, int nz, double * k0 = NULL);
    void getSupercellNestedWannierBands(char * argv0, string fileName, int * size);
  private:
    TBModel * model;
    TBCleanH * ham;
    //BoundaryWilsonH * boundH = NULL;
};

#endif
