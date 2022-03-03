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

void OData::eBands2D(Hamiltonian & ham, int n, int ky, int kz){
  double * k = new double[3];
  k[0] = -M_PI;
  k[1] = ky; 
  k[2] = kz;

  vec ee;
  for(int i = 0; i <= n; i++){
    f << k[0] << " ";
    eig_sym(ee, ham.H(k));
    for(int i = 0; i < size(ee)[0]; i++){
      f << ee.at(i) << " "; 
    }
    f << endl;
    k[0] += 2*M_PI/(double)n;
  }
  delete[] k;
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
