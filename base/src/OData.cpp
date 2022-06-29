#include "OData.h"
#include <iostream>
#include <complex>
#include "Wilson.h"
#include <armadillo>

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

void OData::line(vector<double> line){
  for(int i = 0; i < line.size(); i++){
    f << line[i] << " ";
  }
  if(line.size() != 0){
    f << endl;
  }
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
  //cout << eigVal << endl;
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

void OData::wannierBands(Hamiltonian & ham, int * n, int nWilson, int dirWilson, int nOrbFilled){
  Wilson wilson(&ham);
  wilson.setLoopDir(dirWilson);
  int dim = ham.getNDim();
  double * k = new double[dim];
  for(int i = 0; i < dim; i++){
    k[i] = -M_PI;
  }
  k[dirWilson] = 0.1;
  vec phase(nOrbFilled);


  if(dim == 2){
    double deltaK = 2*M_PI/(double)n[0];
    int xDir; 
    if(dirWilson == 0){
      xDir = 1; 
    }
    else{
      xDir = 0;
    }
    for(int i = 0; i <= n[0]; i++){
      k[xDir] = -M_PI + i*deltaK;
      wilson.setLoopStart(k);
      phase = wilson.wilsonPhases(nWilson, nOrbFilled);
      f << k[xDir] << " " ;
      for(int e = 0; e < size(phase)[0]; e++){
	f << phase[e] << " ";
      }
      f << endl;
    }
  }
  else if(dim == 3){
    double deltaKX = 2*M_PI/(double)n[0];
    double deltaKY = 2*M_PI/(double)n[1];
    int xDir, yDir;
    if(dirWilson == 0){
      xDir = 1;
      yDir = 2;
    }
    else if(dirWilson == 1){
      xDir = 0;
      yDir = 2;
    }
    else{
      xDir = 0; 
      yDir = 1;
    }

    for(int i = 0; i <= n[0]; i++){
      k[xDir] = -M_PI + i*deltaKX;
      for(int e = 0; e <= n[1]; e++){
	k[yDir] = -M_PI + e*deltaKY;
	wilson.setLoopStart(k);
	phase = wilson.wilsonPhases(nWilson, nOrbFilled);
	f << k[xDir] << " " << k[yDir] << " " ;
	for(int j = 0; j < size(phase)[0]; j++){
	  f << phase[j] << " ";
	}
	f << endl;
      }
    }
  }
  else{
    cout << __PRETTY_FUNCTION__ << "\nError: can't produce Wannier bands for a system with " << dim << "dimensions." << endl;
  }

  delete[] k;
}

void OData::nestedWannierBands(Hamiltonian & ham, int nPoints, int xDir, int * nWilson, int * dirWilson, int nOrbFilled){
  Wilson wilson(&ham);
  int dim = ham.getNDim();
  double * k = new double[dim];
  for(int i = 0; i < dim; i++){
    k[i] = 0;
  }
  vec phase;

  if(dim == 3){
    double deltaK = 2*M_PI/(double)nPoints;

    for(int i = 0; i <= nPoints; i++){
      k[xDir] = -M_PI + i*deltaK;
      wilson.setLoopStart(k);
      phase = wilson.nestedWilsonPhases(nWilson, dirWilson, nOrbFilled);
      f << k[xDir] << " " ;
      for(int j = 0; j < size(phase)[0]; j++){
	f << phase[j] << " ";
      }
      f << endl;
    }
  }
  else{
    cout << __PRETTY_FUNCTION__ << "\nError: can't produce nested Wannier bands for a system with " << dim << "dimensions." << endl;
  }

  delete[] k;
}

void OData::supercellWannierBands(Hamiltonian & ham, int * nPoints, int nWilson, int dirWilson, int nOrbFilled){
  Wilson wilson(&ham);
  wilson.setLoopDir(dirWilson);
  int dim = ham.getNDim();
  double * theta = new double[dim];
  for(int i = 0; i < dim; i++){
    theta[i] = ham.getTwists()[i];
  }
  vec phase(nOrbFilled);


  if(dim == 2){
    double deltaTheta = 2*M_PI/(double)nPoints[0];
    int xDir; 
    if(dirWilson == 0){
      xDir = 1; 
    }
    else{
      xDir = 0;
    }
    for(int i = 0; i <= nPoints[0]; i++){
      theta[xDir] = -M_PI + i*deltaTheta;
      ham.setTwists(theta);
      phase = wilson.wilsonPhasesSupercell(nWilson, nOrbFilled);
      f << theta[xDir] << " " ;
      for(int e = 0; e < size(phase)[0]; e++){
	f << phase[e] << " ";
      }
      f << endl;
    }
  }
  else if(dim == 3){
    double deltaThetaX = 2*M_PI/(double)nPoints[0];
    double deltaThetaY = 2*M_PI/(double)nPoints[1];
    int xDir, yDir;
    if(dirWilson == 0){
      xDir = 1;
      yDir = 2;
    }
    else if(dirWilson == 1){
      xDir = 0;
      yDir = 2;
    }
    else{
      xDir = 0; 
      yDir = 1;
    }

    for(int i = 0; i <= nPoints[0]; i++){
      theta[xDir] = -M_PI + i*deltaThetaX;
      for(int e = 0; e <= nPoints[1]; e++){
	theta[yDir] = -M_PI + e*deltaThetaY;
	ham.setTwists(theta);
	phase = wilson.wilsonPhasesSupercell(nWilson, nOrbFilled);
	f << theta[xDir] << " " << theta[yDir] << " " ;
	for(int j = 0; j < size(phase)[0]; j++){
	  f << phase[j] << " ";
	}
	f << endl;
      }
    }
  }
  else{
    cout << __PRETTY_FUNCTION__ << "\nError: can't produce Wannier bands for a system with " << dim << "dimensions." << endl;
  }

  delete[] theta;
}

void OData::supercellNestedWannierBands(Hamiltonian & ham, int nPoints, int xDir, int * nWilson, int * dirWilson, int * nOrbFilled){
  Wilson wilson(&ham);
  int dim = ham.getNDim();
  double * theta = new double[dim];
  for(int i = 0; i < dim; i++){
    theta[i] = ham.getTwists()[i];
  }
  vec phase(nOrbFilled[1]);


  if(dim == 3){
    double deltaTheta = 2*M_PI/(double)nPoints;
    for(int i = 0; i <= nPoints; i++){
      theta[xDir] = -M_PI + i*deltaTheta;
      ham.setTwists(theta);
      phase = wilson.nestedWilsonPhasesSupercell(nWilson, dirWilson, nOrbFilled);
      f << theta[xDir] << " " ;
      for(int e = 0; e < size(phase)[0]; e++){
	f << phase[e] << " ";
      }
      f << endl;
    }
  }
  else{
    cout << __PRETTY_FUNCTION__ << "\nError: can't produce Wannier bands for a system with " << dim << "dimensions." << endl;
  }

  delete[] theta;
}

void OData::matrixWeights(cx_mat & mat){
  double w;
  for(int i = 0; i < size(mat)[0]; i++){
    for(int e = 0; e < size(mat)[1]; e++){
      w = abs(mat(i,e));
      if(w != 0){
	f << e << " " << i << " " << w << endl;
      }
    }
  }
}

void OData::matrixWeightsReal(cx_mat & mat){
  double w;
  for(int i = 0; i < size(mat)[0]; i++){
    for(int e = 0; e < size(mat)[1]; e++){
      w = mat(i,e).real();
      if(w != 0){
	f << e << " " << i << " " << w << endl;
      }
    }
  }
}

void OData::matrixWeightsImag(cx_mat & mat){
  double w;
  for(int i = 0; i < size(mat)[0]; i++){
    for(int e = 0; e < size(mat)[1]; e++){
      w = mat(i,e).imag();
      if(w != 0){
	f << e << " " << i << " " << w << endl;
      }
    }
  }
}
