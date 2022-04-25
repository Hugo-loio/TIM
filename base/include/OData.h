#ifndef ODATA_h
#define ODATA_h

#include <string>
#include <fstream>
#include <vector>
#include "Hamiltonian.h"
using namespace std;

class OData{
  public:
    OData(char * argv0, string fName); //Output argv[0] to argv0
    ~OData();

    void line(string line);
    void data(double ** data, int dim, int nPoints);
    void data(vector<vector<double>> data);

    //x axis is kx, n number of points
    void eBands2D(Hamiltonian & ham, int n, int dir, double * k);
    //x axis is kx, y axis is ky 
    void eBands3D(Hamiltonian & ham, int nx, int ny, int kz = 0);
    //kI defines start and kF defines the finish k point of each segment, segPoints defines number of points in each segment
    void eBandsPath(Hamiltonian & ham, int nSegs, double ** kI, double ** kF, int * segPoints);

    void chargeDensity(Hamiltonian & ham, int nOrb, int nOrbFilled, int * l);

    void wannierBands(Hamiltonian & ham, int * n, int nWilson, int dirWilson, int nOrbFilled); 
    void nestedWannierBands(Hamiltonian & ham, int nPoints, int xDir, int * nWilson, int * dirWilson, int nOrbFilled);

    void supercellWannierBands(Hamiltonian & ham, int * nPoints, int nWilson, int dirWilson, int nOrbFilled);
    void supercellNestedWannierBands(Hamiltonian & ham, int nPoints, int xDir, int * nWilson, int * dirWilson, int * nOrbFilled);

  private:
    ofstream f;
};

#endif
