#include "OData.h"
#include <iostream>

using namespace arma;

OData::OData(char * argv0, string fName){
  string path(argv0);
  string dir = path.substr(0,path.find_last_of('/') + 1);
  f.open(dir + "data/" + fName);
}

OData::~OData(){
  f.close();
}

void OData::line(string line){
  f << line << endl;
}

void OData::data(double ** data, int dim, int nPoints){
  for(int i = 0; i < nPoints; i++){
    for(int e = 0; e < dim; e++){
      f << data[e][i] << " " ; 
    }
    f << endl;
  }
}

void OData::data(vector<vector<double>> data){
  if(data.size() != 0){
    for(int i = 0; i < data[0].size(); i++){
      for(int e = 0; e < data.size(); e++){
	f << data[e][i] << " " ; 
      }
      f << endl;
    }
  }
}

void OData::eBands2D(Hamiltonian & ham, int n, int dir, double * k){
  k[dir] = -M_PI;

  vec ee;
  for(int i = 0; i <= n; i++){
    f << k[dir] << " ";
    eig_sym(ee, ham.H(k));
    for(int i = 0; i < size(ee)[0]; i++){
      f << ee.at(i) << " "; 
    }
    f << endl;
    k[dir] += 2*M_PI/(double)n;
  }
}

void OData::eBands3D(Hamiltonian & ham, int nx, int ny, int kz){
  double * k = new double[3];
  k[0] = -M_PI;
  k[1] = -M_PI; 
  k[2] = kz;

  vec ee;
  for(int e = 0; e <= ny; e++){
    k[0] = -M_PI;
    for(int i = 0; i <= nx; i++){
      f << k[0] << " " << k[1] << " ";
      eig_sym(ee, ham.H(k));
      for(int i = 0; i < size(ee)[0]; i++){
	f << ee.at(i) << " "; 
      }
      f << endl;
      k[0] += 2*M_PI/(double)nx;
    }
    k[1] += 2*M_PI/(double)ny;
  }
  delete[] k;
}

void OData::eBandsPath(Hamiltonian & ham, int nPoints, double ** k, int * segPoints){
  int nDim = ham.getNDim();
  double * kTemp = new double[nDim];

  vec ee;
  for(int i = 0; i < nPoints-1, i++){
    for(int e = 0; e <= segPoints[i]; e++){
      for(int j = 0; j < nDim; i++){
	kTemp[j] = k[i][j]+(double)e*(k[i+1][j]-k[i][j])/(double)segPoints[i];
      }
      f << i+(double)e/segPoints[i] << " ";
      eig_sym(ee, ham.H(kTemp));
      for(int i = 0; i < size(ee)[0]; i++){
	f << ee.at(i) << " "; 
      }
      f << endl;
    }
  }
  delete[] kTemp;
}
