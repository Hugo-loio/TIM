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

void OData::eBandsPath(Hamiltonian & ham, int nSegs, double ** kI, double ** kF, int * segPoints){
  int nDim = ham.getNDim();
  double * kTemp = new double[nDim];

  vec ee;
  for(int i = 0; i < nSegs; i++){
    for(int e = 0; e <= segPoints[i]; e++){
      for(int j = 0; j < nDim; j++){
	kTemp[j] = kI[i][j]+(double)e*(kF[i][j]-kI[i][j])/(double)segPoints[i];
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

void OData::chargeDensity(Hamiltonian & ham, int nOrb, int nOrbFilled, int * l){
  int nDim = ham.getNDim();
  vec eigVal;
  cx_mat eigVec;
  eig_sym(eigVal, eigVec, ham.H());
  int m;
  double rho;

  if(nDim == 1){ 
    cx_mat eigVecOcc = eigVec.cols(0,l[0]*nOrbFilled - 1);
    cx_mat rhoMat = eigVecOcc * (eigVecOcc.t());

    for(int x = 0; x < l[0]; x++){
      rho = 0;
      for(int orb = 0; orb < nOrb; orb++){
	m = orb+x*nOrb; 
	rho += rhoMat(m,m).real();
      }
      f << x + 1 << " " << rho << endl;
    }

  }
  else if(nDim == 2){
    cx_mat eigVecOcc = eigVec.cols(0,l[0]*l[1]*nOrbFilled - 1);
    cx_mat rhoMat = eigVecOcc * (eigVecOcc.t());

    for(int x = 0; x < l[0]; x++){
      for(int y = 0; y < l[1]; y++){
	rho = 0;
	for(int orb = 0; orb < nOrb; orb++){
	  m = orb+x*nOrb + y*nOrb*l[0]; 
	  rho += rhoMat(m,m).real();
	}
	f << x + 1 << " " << y + 1 << " " << rho << endl;
      }
    }
  }
  else if(nDim ==3){
    cx_mat eigVecOcc = eigVec.cols(0,l[0]*l[1]*l[2]*nOrbFilled - 1);
    cx_mat rhoMat = eigVecOcc * (eigVecOcc.t());

    for(int x = 0; x < l[0]; x++){
      for(int y = 0; y < l[1]; y++){
	for(int z = 0; z < l[2]; z++){
	  rho = 0;
	  for(int orb = 0; orb < nOrb; orb++){
	    m = orb+x*nOrb + y*nOrb*l[0] + z*nOrb*l[0]*l[1]; 
	    rho += rhoMat(m,m).real();
	  }
	  f << x + 1 << " " << y + 1 << " " << z + 1 << " " << rho << endl;
	}
      }
    }
  }
}
