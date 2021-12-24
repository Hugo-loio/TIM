#include "TBmod.h"
#include <iostream>

using namespace arma;
using namespace std;

TBmod::TBmod(int ndim, int norb, int * L){
  this->ndim = ndim;
  this->norb = norb;
  this->bc = new int[ndim];
  this->theta = new double[ndim];
  this->L = new int[ndim];
  this->Laccum = new int[ndim];
  this->Lbaccum = new int[ndim];
  this->Lbound = new int[ndim];
  this->nfull_bulk = NULL;
  this->nfull_bound = NULL;
  for(int i = 0; i < ndim; i++){
    this->bc[i] = 1; 
    this->theta[i] = 0;
    this->L[i] = 1;
    this->Lbound[i] = 0;
  }

  set_size(L);
}

TBmod::~TBmod(){
  delete[] L;
  delete[] bc;
  delete[] theta;
  delete[] Laccum;
  delete[] Lbaccum;
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
      if(n[i] > L[i]){
	L[i] = n[i];
      }
    }
  }
  if(newbound){
    set_size(NULL);
  }
}

void TBmod::set_onsite(int norb, complex<double> en){
  os.push_back(Onsite(norb,en));
}

void TBmod::set_bc(int * bc){
  for(int i = 0; i < ndim; i++){
    this->bc[i] = bc[i];
  }
}

void TBmod::set_twists(double * theta){
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
  check_size();
  calc_accum();
  calc_n();
}

void TBmod::calc_accum(){
  //Calculate Laccum
  for(int i = 0; i < ndim; i++){
    Laccum[i] = 1;
    for(int e = 0; e <= i; e++){
      Laccum[i] *= L[e];
    }
  }

  //Calculate Lbaccum
  for(int i = 0; i < ndim; i++){
    Lbaccum[i] = Lbound[0];
    for(int e = 0; e <= i; e++){
      if(Lbound[e] != 0){
	Lbaccum[i] *= Lbound[e];
      }
    }
  }
}

void TBmod::check_size(){
  bool changed = false;
  for(int e = 0; e < hop.size(); e++){
    for(int i = 0; i < ndim; i++){
      if(hop[e].get_n()[i] > L[i]){
	L[i] = hop[e].get_n()[i];
	changed = true;
      }
    }
  }
  if(changed){
    cout << "System size wasn't large enough for the hopping terms\nNew system size:" << endl;
    for(int i = 0; i < ndim; i++){
      cout << "  [" << i << "] = " << L[i] << endl;
    }
  }
}

int TBmod::get_n(int * n){
  int res = 1;
  for(int i = 0; i < ndim; i++){
    if(i == 0){
      res = n[0];
    }
    else{
      res += Laccum[i-1]*n[i];
    }
  }
  return res;
}

void TBmod::calc_n(){
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

void TBmod::next_point_ndim(int depth, int * point, bool up){
  if(depth == ndim-1){

    //increase coordinate
    if(point[ndim-1] == L[ndim-1] - 1){
      point[ndim -1] = 0; 
      next_point_ndim(--depth, point, true);
    }
    else{
      point[ndim-1]++;
    }
  }
  else if(depth != -1){
    if(up){
      //increase coordinate
      if(point[depth] == L[depth] -1 ){
	point[depth] = 0;
	next_point_ndim(--depth, point, true);
      }
      else{
	point[depth]++;
      }
    }
    else{
      //Go to next direction
      next_point_ndim(++depth, point, false);
    }
  }
}

void TBmod::next_point_nrdim(int depth, int * point, bool up){
  if(depth = nrdim-1){
    //increase coordinate
    if(point[nrdim-1] == L[rindex[nrdim -1]] - 1){
      point[nrdim -1] = 0; 
      next_point_nrdim(depth--, point, true);
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
	next_point_nrdim(depth--, point, true);
      }
      else{
	point[depth]++;
      }
    }
    else{
      //Go to next direction
      next_point_nrdim(depth++, point, false);
    }
  }
}

cx_mat TBmod::get_rH(){
  cx_mat res(norb*Laccum[ndim-1], norb*Laccum[ndim-1], fill::zeros);

  complex<double> t;
  int i,j;
  complex<double> phase;
  complex<double> ii(0,1);

  int * nbound = new int[ndim];
  int incn;

  //Hopping terms
  for(int e = 0; e < hop.size(); e++){
    //Increase in n to get neighbour cell in hop (in case of bulk)
    incn = hop[e].get_n()[0];
    for(i = 1; i < ndim; i++){
      incn += hop[e].get_n()[i]*Laccum[i-1];
    }
    //Go through bulk mesh
    for(i = 0; i < Laccum[ndim-1] - Lbaccum[ndim-1]; i++){
      res(hop[e].get_norb1() + norb*nfull_bulk[i][ndim], hop[e].get_norb2() + norb*(nfull_bulk[i][ndim] + incn)) += hop[e].get_hop();
    }
    //Go through boundary mesh and check BCs
    for(i = 0; i < Lbaccum[ndim-1]; i++){
      phase = 1;
      for(j = 0; j < ndim; j++){
	nbound[j] = (nfull_bound[i][j] + hop[e].get_n()[j]) % L[j];
	if(nfull_bound[i][j] + hop[e].get_n()[j] > L[j] - 1){
	  switch(bc[j]){
	    case 0 :
	      phase = 0;
	      break;
	    case 2:
	      phase *= exp(ii*theta[j]);
	  }
	}
      }
      res(hop[e].get_norb1() + norb*nfull_bound[i][ndim], hop[e].get_norb2() + norb*get_n(nbound)) += phase*hop[e].get_hop();
    }
  }

  delete[] nbound;

  res = res + res.t();

  //On-site terms
  int m;
  for(int e  = 0; e < os.size(); e++){
    //Go through bulk mesh
    for(i = 0; i < Laccum[ndim-1] - Lbaccum[ndim-1]; i++){
      m = os[e].get_norb() + norb*nfull_bulk[i][ndim];
      res(m,m) += os[e].get_en();
    }
    //Go through boundary mesh
    for(i = 0; i < Lbaccum[ndim-1]; i++){
      m = os[e].get_norb() + norb*nfull_bound[i][ndim];
      res(m,m) += os[e].get_en();
    }
  }

  return res;
}
