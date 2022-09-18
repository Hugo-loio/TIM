#ifndef ENGAP_H
#define ENGAP_H

#include "Hamiltonian.h"

using namespace arma;
using namespace std;

class EnGap{
  public:
    EnGap(Hamiltonian * ham);
    ~EnGap();

    double getGap(double en = 0, double * k = NULL);

  private:
    Hamiltonian * ham;
};

#endif
