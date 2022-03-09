#ifndef BBH_H
#define BBH_H

#include "TBCleanH.h"
#include <armadillo>
#include <string>
#include "Wilson.h"
#include "BoundaryWilsonH.h"

using namespace arma;
using namespace std;

class BBH2D{
  public:
    BBH2D(double t1 = 1, double t2 = 2);
    ~BBH2D();

    void setIntraHop(double);
    void setInterHop(double);
    //void setSparse(bool);

    double berryPhase(int n, int dir = 0, double * k0 = NULL);
    void getBands(char * argv0, string fileName, int nx, int ny);
    void getWannierBands(char * argv0, string fileName, int dir, int n);
    //double getQuadrupoleNested(int nx, int ny, double * k0 = NULL);
  private:
    TBModel * model;
    TBCleanH * ham;
    BoundaryWilsonH * boundH;
};

#endif
