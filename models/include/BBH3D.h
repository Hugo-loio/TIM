#ifndef BBH3D_H
#define BBH3D_H

#include "TBCleanH.h"
#include <armadillo>
#include <string>
#include "Wilson.h"
#include "BoundaryWilsonH.h"

using namespace arma;
using namespace std;

class BBH3D{
  public:
    //Intracell hopping t1, intercell hopping t2
    BBH3D(double t1 = 1, double t2 = 2);
    ~BBH3D();

    void setIntraHop(double);
    void setInterHop(double);
    //void setSparse(bool);

    //double berryPhase(int n, int dir = 0, double * k0 = NULL);
    void getBands(char * argv0, string fileName);
    //void getWannierBands(char * argv0, string fileName, int dir, int n);
    cx_mat getH(double * k = NULL);
    //double getQuadrupoleNested(int nx, int ny, double * k0 = NULL);
  private:
    TBModel * model;
    TBCleanH * ham;
    //BoundaryWilsonH * boundH = NULL;
};

#endif
