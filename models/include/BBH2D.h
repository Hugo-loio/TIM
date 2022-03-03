#ifndef BBH_H
#define BBH_H

#include "TBCleanH.h"
#include <armadillo>
#include <string>
#include "Wilson.h"

using namespace arma;
using namespace std;

class BBH2D{
  public:
    BBH2D(double t1 = 1, double t2 = 2);
    ~BBH2D();

    void setIntraHop(double);
    void setInterHop(double);
    void setSparse(bool);

    double berryPhase(int n, int dir = 0, double * k0 = NULL);
    void getBands(char * arv0, string fileName, int nx, int ny);
  private:
    TBModel * model;
    TBCleanH * ham;
};

#endif
