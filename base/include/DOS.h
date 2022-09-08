#ifndef DOS_H
#define DOS_H

#include "Hamiltonian.h"

using namespace std;

class DOS{
  public:
    DOS(Hamiltonian * ham);
    ~DOS();

    double kpm(double en, int nMoments, int nRandVecs, double * k = NULL);

  private:
    Hamiltonian * ham;
    double a,b;
    double err = 0.01;
    double * mu;
    int nMoments;
    int nRandVecs;
    bool momentsFound = false;
    bool rescalingFound = false;

    void findRescaling(double * k = NULL);
    void calculateMoments(double * k = NULL);
    void randomize(cx_vec & rand, int d);
    double jacksonKernel(int n);
    double chebyshev(int n, double x);
};

#endif
