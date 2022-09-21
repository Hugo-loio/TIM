#ifndef LOCALIZATIONPROPS_H
#define LOCALIZATIONPROPS_H

#include "Hamiltonian.h"

using namespace std;

class LocalizationProps{
  public:
    LocalizationProps(Hamiltonian * ham);
    ~LocalizationProps();

    //Inverse Participation Ratio
    double ipr(int vol, int nOrb, int nStates, double en, double * k = NULL);
    //Level spacing statistics
    double lsr(int nStates, double en, double * k = NULL);
    //Energy gap
    double gap(double en = 0, double * k = NULL);
    //Transfer Matrix Method, localization length
    double tmm(int nLayers, int qrIt, double en, double * k = NULL);

    void setForceDiag(){forceDiag = true;}

  private:
    Hamiltonian * ham;
    int maxItTMM = 1E6;
    double tmmErr = 0.01;
    cx_mat eigVec;
    vec eigVal;
    int nStatesFound = 0;
    double lastEn;
    bool forceDiag = false;

    void sparseDiag(int nStates, double en, double * k);
    bool testTmmConv(double * c, double * d, int size, int nQR, int & minIndex);
};

#endif
