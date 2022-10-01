#ifndef DOS_H
#define DOS_H

#include "Hamiltonian.h"

using namespace std;

class DOS{
  public:
    DOS(Hamiltonian * ham);
    ~DOS();

    double kpm(double en, int nMoments, int nRandVecs, double * k = NULL);
    void setKpmERange(double eMin, double eMax);

  private:
    Hamiltonian * ham;
    double a,b;
    double err = 0.01;
    double * mu;
    int nMoments;
    int nRandVecs;
    bool momentsFound = false;
    bool rescalingFound = false;
    bool customERange = false;

    void findRescaling(double * k = NULL);
    template <class mat> void calculateMoments(mat h);
    void randomize(cx_vec & rand, int d);
    double jacksonKernel(int n);
    double chebyshev(int n, double x);
};

#endif
