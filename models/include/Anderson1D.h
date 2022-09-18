#ifndef ANDERSON1D_H
#define ANDERSON1D_H

#include "DisorderedAnderson1D.h"

using namespace std;

class Anderson1D{
  public:
    Anderson1D(double t = -1);
    ~Anderson1D();

    void setHop(double);
    void setSize(int * l);
    void setW(double);

    void generateDisorder();

    double ipr(int nStates, double en = 0);
    double getTMM(int qrIt, double en);

    cx_mat getHam();
    void test(char * argv0 = NULL);

  private:
    TBModel * model;
    DisorderedAnderson1D * ham;
};

#endif
