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

    void findRescaling(double * k = NULL);
    void calculateMoments(double * k = NULL);
    double jacksonKernel(int n);
};

#endif
