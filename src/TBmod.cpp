#include "TBmod.h"
#include <iostream>

using namespace arma;
using namespace std;

TBmod::TBmod(int ndim, int norb, int * L){
  this->ndim = ndim;
  this->norb = norb;
  this->pcb = new bool[ndim];
  this->theta = new double[ndim];
  this->L = new int[ndim];
  this->Laccum = new int[ndim];
  this->Lbound = new int[ndim];
  for(int i = 0; i < ndim; i++){
    this->pcb[i] = true; 
    this->theta[i] = 0;
    this->L[i] = 1;
    this->Lbound = 0;
  }
  set_size(L);
}

TBmod::~TBmod(){
  delete[] L;
  delete[] pcb;
  delete[] theta;
  delete[] Laccum;
  delete[] Lbound;
  delete[] nfull_bulk;
  delete[] nfull_bound;
  //delete[] n;
}

void TBmod::set_hop(int norb1, int norb2, int * n, complex<double> hop){
  this->hop.push_back(Hop(norb1, norb2, n, hop, ndim));
  bool newbound = false;
  for(int i = 0; i < ndim; i++){
    if(n[i] > Lbound[i]){
      newbound = true;
      Lbound[i] = n[i];
    }
  }
  if(newbound){
    calc_n();
  }
}

void TBmod::set_onsite(int norb, complex<double> en){
  os.push_back(Onsite(norb,en));
}

void TBmod::set_periodic(bool * pbc){
  for(int i = 0; i < ndim; i++){
    this->pcb[i] = pbc[i];
  }
}

void TBmod::set_twisted(double * theta){
  for(int i = 0; i < ndim; i++){
    this->theta[i] = theta[i];
  }
}

void TBmod::set_size(int * L){
  if(L != NULL){
    for(int i = 0; i < ndim; i++){
      this->L[i] = L[i];
    }
  }

  //Calculate Laccum
  for(int i = 0; i < ndim; i++){
    Laccum[i] = 1;
    for(int e = 0; e <= i; e++){
      Laccum[i] *= L[e];
    }
  }

  //calculate system volume
  vol = 1;
  for(int i = 0; i < ndim; i++){
    vol *= L[i];
  }

  calc_n();
}


bool TBmod::check_size(){
  bool res = true;
  for(int e = 0; e < hop.size(); e++){
    for(int i = 0; i < ndim; i++){
      if(!pcb[i]){
	if(hop[e].get_n()[i] > L[i]){
	  res = false;
	}
      }
    }
  }
  return res;
}

void TBmod::calc_n(){
  delete[] nfull_bulk;
  nfull_bulk = new int*[Laccum[ndim-1]];
  /*
     int e,k;
     for(int i = 0; i < Laccum[ndim-1]; i++){
     nfull[i] = new int[ndim];
     for(e = 0; e < ndim; e++){
     nfull[i][e] = i;
     for(k = 0; k < e; k++){
     nfull[i][e] -= nfull[i][k];
     }
     nfull[i][e] = nfull[i][e] % Laccum[e];
     if(e > 0){
     nfull[i][e] = nfull[i][e]/Laccum[e-1];
     }
     }
     }
     */
}

void TBmod::next_point_full(int i, int depth, int * point, bool up){
  if(depth = ndim-1){
    /*
    //TODO:: check if boundary or bulk
    //fill nfull
    int n = 1;
    for(int e = 0; e < ndim; e++){
    if(e == 0){
    n = val[0];
    }
    else{
    n += Laccum[e-1]*val[e];
    }
    nfull_bulk[i][e] = val[e];
    }
    nfull_bulk[ndim] = n;
    */

    //increase coordinate
    if(point[ndim-1] == L[ndim-1] - 1){
      point[ndim -1] = 0; 
      loop_nfull(i, depth--, point, true);
    }
    else{
      point[ndim-1]++;
    }
  }
  else{
    if(up){
      //increase coordinate
      if(point[depth] == L[depth] -1 ){
	point[depth] = 0;
	loop_nfull(i, depth--, point, true);
      }
      else{
	point[depth]++;
      }
    }
    else{
      loop_nfull(i, depth++, point, false);
    }
  }
}

/*
   int TBmod::get_m1(int * n, Hop & hop){
   int * temp = new int[ndim];
   int temp2 = 0;

   int e = 0;
   for(int i = 0; i < ndim; i++){
   temp[i] = norb*n[i];
   for(e = 0; e < i; e++){
   temp[i] *= L[e];
   }
   temp2 += temp[i];
   }
   delete[] temp;

   return hop.get_norb1() + temp2;
   }

   int TBmod::get_m2(int * n, Hop & hop){
   int * temp = new int[ndim];
   int temp2 = 0;

   int e = 0;
   for(int i = 0; i < ndim; i++){
   temp[i] = norb*((n[i] + hop.get_n()[i]) % L[i]);
   for(e = 0; e < i; e++){
   temp[i] *= L[e];
   }
   temp2 += temp[i];
   }
   delete[] temp;

   return hop.get_norb2() + temp2;
   }
   */

cx_mat TBmod::get_rH(){
  cx_mat res(norb*Laccum[ndim-1], norb*Laccum[ndim-1], fill::zeros);

  if(check_size()){
    complex<double> t;

    int i,j;
    complex<double> phase;
    complex<double> ii(0,1);
    for(int e = 0; e < hop.size(); e++){
      t = hop[e].get_hop();
      for(i = 0; i < Laccum[ndim -1] ; i++){
	phase = 1;
	for(j = 0; j < ndim; j++){
	  if(nfull[i][j] + hop[e].get_n()[j] > L[j] - 1){
	    if(theta[j] == 0 && !pcb[j]){
	      phase = 0;
	    }
	    else{
	      phase *= exp(ii*theta[j]);
	    }
	  }
	}
	res(get_m1(nfull[i],hop[e]),get_m2(nfull[i],hop[e])) += t*phase;
      }
    }
  }

  else{
    cout << "Error: Hopping terms aren't compatible with system size" << endl;
  }
  return res + res.t();
}
