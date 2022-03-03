#ifndef SSH2D_H
#define SSH2D_H

#include "TBCleanH.h"
#include <armadillo>
#include <string>
#include "Wilson.h"

using namespace arma;
using namespace std;

class SSH2D{
  public:
    SSH2D(double t1 = 1, double t2 = 2);
    ~SSH2D();

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
