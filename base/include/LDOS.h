#ifndef LDOS_H
#define LDOS_H

#include "Hamiltonian.h"

using namespace std;

class LDOS{
  public:
    LDOS(Hamiltonian * ham, int startIndex, int endIndex);
    ~LDOS();

    double kpm(double en, int nMoments, double * k = NULL);
    void setKpmERange(double eMin, double eMax);

    double diag(double en, double range, double * k = NULL);

  private:
    Hamiltonian * ham;
    //Unit cell vector
    int startIndex;
    int endIndex;
    double a,b;
    double err = 0.01;
    double * mu;
    int nMoments;
    bool momentsFound = false;
    bool rescalingFound = false;
    bool customERange = false;

    void findRescaling(double * k = NULL);
    template <class mat> void calculateMoments(mat h);
    double jacksonKernel(int n);
    double chebyshev(int n, double x);
};

#endif
