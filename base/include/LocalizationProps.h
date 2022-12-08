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
    //Transfer Matrix Method, localization length + relative error
    vector<double> tmm(int nLayers, int qrIt, double en, double * k = NULL);
    vector<double> tmmSpecial(int nLayers, int qrIt, double en, double * k = NULL);
    vector<double> tmmSpecial2(int nLayers, int qrIt, double en, double * k = NULL);
    vector<double> tmmSpecialReal(int nLayers, int qrIt, double en, double * k = NULL);
    vector<double> tmmSpecialReal2(int nLayers, int qrIt, double en, double * k = NULL);
    vector<double> tmmSpecialReal3(int nLayers, int qrIt, double en, double * k = NULL);

    void setForceDiag(){forceDiag = true;}
  private:
    Hamiltonian * ham;
    int maxItTMM = 1E7;
    double tmmErr = 0.02;
    cx_mat eigVec;
    vec eigVal;
    int nStatesFound = 0;
    double lastEn;
    bool forceDiag = false;

    void sparseDiag(int nStates, double en, double * k);
    // Auxiliar to TMM
    bool testTmmConv(double * c, double * d, int size, int nQR, int & minIndex, double & err);
    template <class mat> void diagInverse(mat & v, int size);
    template <class mat> mat diagLeftMult(mat & d, mat m, int size);
    template <class mat> mat diagDoubleMult(mat & d1, mat d2, int size);
    template <class mat> mat updateT(mat & t, mat & tAux, int size);
};

#endif
