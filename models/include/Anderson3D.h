#ifndef ANDERSON3D_H
#define ANDERSON3D_H

#include "DisorderedAnderson3D.h"

using namespace std;

class Anderson3D{
  public:
    Anderson3D(double t = -1);
    ~Anderson3D();

    void setHop(double);
    void setSize(int * l);
    void setW(double);

    void generateDisorder();

    double ipr(int nStates, double en = 0);
    vector<double> getTMM(int qrIt, double en, int m);

    cx_mat getHam();
    void test(char * argv0 = NULL);

  private:
    TBModel * model;
    DisorderedAnderson3D * ham;
};

#endif
