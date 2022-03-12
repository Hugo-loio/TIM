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
    //nPoints dimensional k defines path, (nPoints-1) dimensional segPoints defines number of points in each segment
    void eBandsPath(Hamiltonian & ham, int nPoints, double ** k, int * segPoints);

  private:
    ofstream f;
};

#endif
