#ifndef SSH_H
#define SSH_H

#include "TBCleanH.h"
#include <armadillo>
#include <string>
#include "Wilson.h"

using namespace arma;
using namespace std;

class SSH{
  public:
    SSH(double t1 = 1, double t2 = 2);
    ~SSH();

    void setIntraHop(double);
    void setInterHop(double);

    double berryPhase(int n);
    void getBands(char * arv0, string fileName, int n);
  private:
    TBModel * model;
    TBCleanH * ham;
};

#endif
