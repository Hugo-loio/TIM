#ifndef SOTAI_H
#define SOTAI_H

#include "TBCleanH.h"
#include <armadillo>
#include <string>
#include "Wilson.h"

using namespace arma;
using namespace std;

class SOTAI{
  public:
    //Intracell hopping t1, intercell hopping t2
    SOTAI(double m);
    ~SOTAI();

    void setM(double);

    double berryPhase(int n, int dir = 0, double * k0 = NULL);
    void getBands(char * argv0, string fileName, int nx, int ny);
    void getWannierBands(char * argv0, string fileName, int dir);
    void getChargeDensity(char * argv0, string fileName, int nx, int ny, int nOrbFilled);
    //cx_mat getH(double * k = NULL);
    double getQuadrupoleNested(int nx, int ny, double * k0 = NULL);
    double getQuadrupoleNestedSupercell(int * l, int * n);

    void getSupercellWannierBands(char * argv0, string fileName, int nx, int ny, int dirWilson);

  private:
    TBModel * model;
    TBCleanH * ham;
};

#endif
