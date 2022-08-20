#ifndef LOCALIZATIONSTATS_H
#define LOCALIZATIONSTATS_H

#include "Hamiltonian.h"

using namespace std;

class LocalizationStats{
  public:
    LocalizationStats(Hamiltonian * ham);
    ~LocalizationStats();

    //Inverse Participation Ratio
    double ipr(int vol, int nOrb, int nStates, double * k = NULL);
    //Transfer Matrix Method, localization length
    double tmm(int nLayers, int qrIt, double en, double * k = NULL);
  private:
    Hamiltonian * ham;
    int maxItTMM = 1E6;
    double tmmErr = 0.01;

    bool testTmmConv(double * c, double * d, int size, int nQR, int & minIndex);
};

#endif
