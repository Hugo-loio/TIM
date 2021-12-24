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

  //Calculate Lbacuum
  for(int i = 0; i < ndim; i++){
    Lbaccum[i] = 1;
    for(int e = 0; e <= i; e++){
      Lbaccum[i] *= L[e]-Lbound[e];
    }
  }

  //Calculate 

  /*
  //calculate system volume
  vol = 1;
  for(int i = 0; i < ndim; i++){
  vol *= L[i];
  }
  */

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

int TBmod::get_n(int * n){
  int res = 1;
  for(int i = 0; i < ndim; i++){
    if(i == 0){
      res = n[0];
    }
    else{
      res += Laccum[e-1]*n[e];
    }
  }
  return res;
}

void TBmod::calc_n(){
  if(check_size()){
    if(nfull_bulk != NULL){
      delete[] nfull_bulk;
    }
    if(nfull_bound != NULL){
      delete[] nfull_bound;
    }

    nfull_bulk = new int * [Laccum[ndim-1]-Lbaccum[ndim-1]];
    nfull_bound = new int * [Lbaccum[ndim -1]];

    int * point = new int[ndim];
    for(int i = 0; i < ndim; i++){
      point[i] = 0;
    }

    int e,j;
    //nfull
    int count_bulk = 0;
    int count_bound = 0;
    bool bound = false;
    for(int i = 0; i < Laccum[ndim-1];i++){
      for(e = 0; e < ndim; e++){
	if(point[e] > L[e] - Lbound[e] - 1){
	  bound = true;
	}
      }
      if(bound){
	nfull_bound[count_bound] = new int[ndim + 1];
	for(j = 0; j < ndim; j++){
	  nfull_bound[count_bound][j] = point[j];
	}
	nfull_bound[count_bound][ndim] = get_n(point);
	count_bound++;
      }
      else{
	nfull_bulk[count_bulk] = new int[ndim + 1];
	for(j = 0; j < ndim; j++){
	  nfull_bulk[count_bulk][j] = point[j];
	}
	nfull_bulk[count_bulk][ndim] = get_n(point);
	count_bulk++;
      }
      bound = false;
      next_point_ndim(0, point, false);
    }
  }
  else{
    cout << "Error: Hopping terms aren't compatible with system size" << endl;
  }
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

void TBmod::next_point_ndim(int depth, int * point, bool up){
  if(depth = ndim-1){

    //increase coordinate
    if(point[ndim-1] == L[ndim-1] - 1){
      point[ndim -1] = 0; 
      next_point_ndim(depth--, point, true);
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
	next_point_ndim(depth--, point, true);
      }
      else{
	point[depth]++;
      }
    }
    else{
      //Go to next direction
      next_point_ndim(depth++, point, false);
    }
  }
}

void TBmod::next_point_nrdim(int i, int depth, int * point, bool up){
  if(depth = nrdim-1){
    //increase coordinate
    if(point[nrdim-1] == L[rindex[nrdim -1]] - 1){
      point[nrdim -1] = 0; 
      next_point_nrdim(i, depth--, point, true);
    }
    else{
      point[nrdim-1]++;
    }
  }
  else{
    if(up){
      //increase coordinate
      if(point[depth] == L[rindex[depth]] -1 ){
	point[depth] = 0;
	next_point_nrdim(i, depth--, point, true);
      }
      else{
	point[depth]++;
      }
    }
    else{
      //Go to next direction
      next_point_ndim(i, depth++, point, false);
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
