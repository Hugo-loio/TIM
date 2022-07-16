#ifndef LOCALIZATIONSTATS_H
#define LOCALIZATIONSTATS_H

#include "Hamiltonian.h"

using namespace std;

class LocalizationStats{
  public:
    LocalizationStats(Hamiltonian * ham);
    ~LocalizationStats();

    //double ipr();

  private:
    Hamiltonian * ham;
};

#endif
