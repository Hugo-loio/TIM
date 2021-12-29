#include "TBmod.h"
#include <iostream>

using namespace arma;
using namespace std;

TBmod::TBmod(int ndim, int norb, int * L){
  this->ndim = ndim;
  this->nrdim = 0;
  this->norb = norb;
  this->bc = new int[ndim];
  this->theta = new double[ndim];
  this->L = new int[ndim];
  this->Laccum = new int[ndim];
  this->Lbound = new int * [ndim];
  this->nfull_bulk = NULL;
  this->nfull_bound = NULL;
  this->n_bulk = NULL;
  this->n_bound = NULL;
  this->rindex = new int[nrdim];
  this->Lraccum = new int[nrdim];
  for(int i = 0; i < ndim; i++){
    this->bc[i] = 1; 
    this->theta[i] = 0;
    this->L[i] = 1;
    this->Lbound[i] = new int[2];
    this->Lbound[i][0] = 0;
    this->Lbound[i][1] = 0;
  }

  set_size(L);
}

TBmod::~TBmod(){
  delete[] L;
  delete[] bc;
  delete[] theta;
  delete[] Lraccum;
  delete[] Laccum;
  delete[] Lbound;
  delete[] nfull_bulk;
  delete[] nfull_bound;
  delete[] n_bulk;
  delete[] n_bound;
  delete[] rindex;
}

void TBmod::set_hop(int norb1, int norb2, int * n, complex<double> hop){
  this->hop.push_back(Hop(norb1, norb2, n, hop, ndim));
  bool newbound = false;
  for(int i = 0; i < ndim; i++){
    if(abs(n[i]) > L[i]){
      L[i] = n[i];
    }
    if(n[i] > 0 && n[i] > Lbound[i][1]){
      newbound = true;
      Lbound[i][1] = n[i];
    }
    else if(n[i] < 0 && n[i] < -Lbound[i][0]){
      newbound = true;
      Lbound[i][0] = -n[i];
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
  nrdim = 0;
  for(int i = 0; i < ndim; i++){
    this->bc[i] = bc[i];
    if(bc[i] != 1){
      nrdim++;
    }
  }

  delete[] rindex;
  rindex = new int [nrdim];
  int e = 0;
  for(int i = 0; i < ndim; i++){
    if(bc[i] != 1){
      rindex[e] = i;
      e++;
    }
  }

  delete[] Lraccum;
  Lraccum = new int [nrdim];

  calc_vol();
  calc_n();
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
  calc_vol();
  calc_nfull();
  calc_n();
}

void TBmod::calc_vol(){
  //Calculate Laccum
  for(int i = 0; i < ndim; i++){
    Laccum[i] = 1;
    for(int e = 0; e <= i; e++){
      Laccum[i] *= L[e];
    }
  }

  //Calculate Lraccum
  for(int i = 0; i < nrdim; i++){
    Lraccum[i] = 1;
    for(int e = 0; e <= i; e++){
      Lraccum[i] *= L[rindex[e]];
    }
  }

  //Calculate Vbound and Vbulk for full mesh
  Vbulkfull = 1;
  Vboundfull = 0;
  for(int i = 0; i < ndim; i++){
    Vbulkfull *= L[i] - (Lbound[i][0] + Lbound[i][1]);
  }
  if(Vbulkfull < 0){
    Vbulkfull = 0;
  }
  Vboundfull = Laccum[ndim -1] - Vbulkfull;

  //Calculate Vbound and Vbulk for reduced mesh
  Vbulk = 1;
  Vbound = 0;
  int e;
  for(int i = 0; i < nrdim; i++){
    e = rindex[i];
    Vbulk *= L[e] - (Lbound[e][0] + Lbound[e][1]);
  }
  if(Vbulk <  0){
    Vbulk = 0;
  }
  if(nrdim != 0){
    Vbound = Lraccum[nrdim -1] - Vbulk;
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

int TBmod::get_nfull(int * n){
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

int TBmod::get_n(int * n){
  int res = 1;
  for(int i = 0; i < nrdim; i++){
    if(i == 0){
      res = n[0];
    }
    else{
      res += Lraccum[i-1]*n[i];
    }
  }
  return res;
}

void TBmod::calc_nfull(){
  if(nfull_bulk != NULL){
    delete[] nfull_bulk;
  }
  if(nfull_bound != NULL){
    delete[] nfull_bound;
  }

  nfull_bulk = new int * [Vbulkfull];
  nfull_bound = new int * [Vboundfull];

  int * point = new int[ndim];
  for(int i = 0; i < ndim; i++){
    point[i] = 0;
  }

  int e,j;
  int count_bulk = 0;
  int count_bound = 0;
  bool bound = false;
  for(int i = 0; i < Laccum[ndim-1];i++){
    for(e = 0; e < ndim; e++){
      if(point[e] > L[e] - Lbound[e][1] - 1 || point[e] < Lbound[e][0]){
	bound = true;
      }
    }
    if(bound){
      nfull_bound[count_bound] = new int[ndim + 1];
      for(j = 0; j < ndim; j++){
	nfull_bound[count_bound][j] = point[j];
      }
      nfull_bound[count_bound][ndim] = get_nfull(point);
      count_bound++;
    }
    else{
      nfull_bulk[count_bulk] = new int[ndim + 1];
      for(j = 0; j < ndim; j++){
	nfull_bulk[count_bulk][j] = point[j];
      }
      nfull_bulk[count_bulk][ndim] = get_nfull(point);
      count_bulk++;
    }
    bound = false;
    next_point_nfull(0, point, false);
  }
}

void TBmod::calc_n(){
  if(n_bulk != NULL){
    delete[] n_bulk;
  }
  if(n_bound != NULL){
    delete[] n_bound;
  }

  n_bulk = new int * [Vbulk];
  n_bound = new int * [Vbound];

  int * point = new int[nrdim];
  for(int i = 0; i < nrdim; i++){
    point[i] = 0;
  }

  int e,j;
  int count_bulk = 0;
  int count_bound = 0;
  bool bound = false;
  if(nrdim != 0){
    for(int i = 0; i < Lraccum[nrdim -1] ;i++){
      for(e = 0; e < nrdim; e++){
	j = rindex[e];
	if(point[e] > L[j] - Lbound[j][1] - 1 || point[e] < Lbound[j][0]){
	  bound = true;
	}
      }
      if(bound){
	n_bound[count_bound] = new int[nrdim + 1];
	for(e = 0; e < nrdim; e++){
	  n_bound[count_bound][e] = point[e];
	}
	n_bound[count_bound][nrdim] = get_n(point);
	count_bound++;
      }
      else{
	n_bulk[count_bulk] = new int[nrdim + 1];
	for(e = 0; e < nrdim; e++){
	  n_bulk[count_bulk][e] = point[e];
	}
	n_bulk[count_bulk][nrdim] = get_n(point);
	count_bulk++;
      }
      bound = false;
      next_point_n(0, point, false);
    }
  }
  else{
    n_bulk[0] = new int[2];
    n_bulk[0][0] = 0;
    n_bulk[0][1] = 0;
  }
}

void TBmod::next_point_nfull(int depth, int * point, bool up){
  if(depth == ndim-1){

    //increase coordinate
    if(point[ndim-1] == L[ndim-1] - 1){
      point[ndim -1] = 0; 
      next_point_nfull(--depth, point, true);
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
	next_point_nfull(--depth, point, true);
      }
      else{
	point[depth]++;
      }
    }
    else{
      //Go to next direction
      next_point_nfull(++depth, point, false);
    }
  }
}

void TBmod::next_point_n(int depth, int * point, bool up){
  if(depth == nrdim-1){
    //increase coordinate
    if(point[nrdim-1] == L[rindex[nrdim -1]] - 1){
      point[nrdim -1] = 0; 
      next_point_n(--depth, point, true);
    }
    else{
      point[nrdim-1]++;
    }
  }
  else if(depth != -1){
    if(up){
      //increase coordinate
      if(point[depth] == L[rindex[depth]] -1 ){
	point[depth] = 0;
	next_point_n(--depth, point, true);
      }
      else{
	point[depth]++;
      }
    }
    else{
      //Go to next direction
      next_point_n(++depth, point, false);
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
    for(i = 0; i < Vbulkfull; i++){
      res(hop[e].get_norb1() + norb*nfull_bulk[i][ndim], hop[e].get_norb2() + norb*(nfull_bulk[i][ndim] + incn)) += hop[e].get_hop();
    }
    //Go through boundary mesh and check BCs
    for(i = 0; i < Vboundfull; i++){
      phase = 1;
      for(j = 0; j < ndim; j++){
	nbound[j] = (nfull_bound[i][j] + hop[e].get_n()[j] + L[j]) % L[j];
	if(nfull_bound[i][j] + hop[e].get_n()[j] > L[j] - 1){
	  switch(bc[j]){
	    case 0 :
	      phase = 0;
	      break;
	    case 2:
	      phase *= exp(-ii*theta[j]);
	  }
	}
	else if(nfull_bound[i][j] + hop[e].get_n()[j] < 0){
	  switch(bc[j]){
	    case 0 :
	      phase = 0;
	      break;
	    case 2:
	      phase *= exp(ii*theta[j]);
	  }
	}
      }
      res(hop[e].get_norb1() + norb*nfull_bound[i][ndim], hop[e].get_norb2() + norb*get_nfull(nbound)) += phase*hop[e].get_hop();
    }
  }

  delete[] nbound;

  res = res + res.t();

  //On-site terms
  int m;
  for(int e  = 0; e < os.size(); e++){
    //Go through bulk mesh
    for(i = 0; i < Vbulkfull; i++){
      m = os[e].get_norb() + norb*nfull_bulk[i][ndim];
      res(m,m) += os[e].get_en();
    }
    //Go through boundary mesh
    for(i = 0; i < Vboundfull; i++){
      m = os[e].get_norb() + norb*nfull_bound[i][ndim];
      res(m,m) += os[e].get_en();
    }
  }

  return res;
}

sp_cx_mat TBmod::get_sprH(){
  sp_cx_mat res(norb*Laccum[ndim-1], norb*Laccum[ndim-1]);

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
    for(i = 0; i < Vbulkfull; i++){
      res(hop[e].get_norb1() + norb*nfull_bulk[i][ndim], hop[e].get_norb2() + norb*(nfull_bulk[i][ndim] + incn)) += hop[e].get_hop();
    }
    //Go through boundary mesh and check BCs
    for(i = 0; i < Vboundfull; i++){
      phase = 1;
      for(j = 0; j < ndim; j++){
	nbound[j] = (nfull_bound[i][j] + hop[e].get_n()[j] + L[j]) % L[j];
	if(nfull_bound[i][j] + hop[e].get_n()[j] > L[j] - 1){
	  switch(bc[j]){
	    case 0 :
	      phase = 0;
	      break;
	    case 2:
	      phase *= exp(-ii*theta[j]);
	  }
	}
	else if(nfull_bound[i][j] + hop[e].get_n()[j] < 0){
	  switch(bc[j]){
	    case 0 :
	      phase = 0;
	      break;
	    case 2:
	      phase *= exp(ii*theta[j]);
	  }
	}
      }
      res(hop[e].get_norb1() + norb*nfull_bound[i][ndim], hop[e].get_norb2() + norb*get_nfull(nbound)) += phase*hop[e].get_hop();
    }
  }

  delete[] nbound;

  res = res + res.t();

  //On-site terms
  int m;
  for(int e  = 0; e < os.size(); e++){
    //Go through bulk mesh
    for(i = 0; i < Vbulkfull; i++){
      m = os[e].get_norb() + norb*nfull_bulk[i][ndim];
      res(m,m) += os[e].get_en();
    }
    //Go through boundary mesh
    for(i = 0; i < Vboundfull; i++){
      m = os[e].get_norb() + norb*nfull_bound[i][ndim];
      res(m,m) += os[e].get_en();
    }
  }

  return res;
}

cx_mat TBmod::get_H(double * k){
  int size; 

  if(nrdim != 0){
    size = norb*Lraccum[nrdim-1];
  }
  else{
    size = norb;
  }
  cx_mat res(size, size, fill::zeros);

  complex<double> t;
  int i,j,l;
  complex<double> phase;
  complex<double> kphase;
  complex<double> ii(0,1);

  int * nbound = new int[nrdim];
  int incn;

  //Hopping terms
  for(int e = 0; e < hop.size(); e++){
    //Increase in n to get neighbour cell in hop (in case of bulk)
    incn = 0;
    j = 0;
    kphase = 1;
    for(i = 0; i < ndim; i++){
      if(bc[i] != 1){
	if(i == rindex[0]){
	  incn = hop[e].get_n()[i];
	}
	else{
	  incn += hop[e].get_n()[i]*Lraccum[j-1];
	}
	j++;
      }
      else{
	//Extra phase due to k space hamiltonian
	kphase *= exp(ii*k[i]*(double)hop[e].get_n()[i]);
      }
    }

    //Go through bulk mesh
    for(i = 0; i < Vbulk; i++){
      res(hop[e].get_norb1() + norb*n_bulk[i][nrdim], hop[e].get_norb2() + norb*(n_bulk[i][nrdim] + incn)) += hop[e].get_hop()*kphase;
    }
    //Go through boundary mesh and check BCs
    for(i = 0; i < Vbound; i++){
      phase = 1;
      for(j = 0; j < nrdim; j++){
	l = rindex[j];
	nbound[j] = (n_bound[i][j] + hop[e].get_n()[l] + L[l]) % L[l];
	if(n_bound[i][j] + hop[e].get_n()[l] > L[l] - 1){
	  switch(bc[l]){
	    case 0 :
	      phase = 0;
	      break;
	    case 2:
	      phase *= exp(-ii*theta[j]);
	  }
	}
	else if(nfull_bound[i][j] + hop[e].get_n()[l] < 0){
	  switch(bc[j]){
	    case 0 :
	      phase = 0;
	      break;
	    case 2:
	      phase *= exp(ii*theta[j]);
	  }
	}
      }
      res(hop[e].get_norb1() + norb*n_bound[i][nrdim], hop[e].get_norb2() + norb*get_n(nbound)) += kphase*phase*hop[e].get_hop();
    }
  }

  delete[] nbound;

  res = res + res.t();

  //On-site terms
  int m;
  for(int e  = 0; e < os.size(); e++){
    //Go through bulk mesh
    for(i = 0; i < Vbulk; i++){
      m = os[e].get_norb() + norb*n_bulk[i][nrdim];
      res(m,m) += os[e].get_en();
    }
    //Go through boundary mesh
    for(i = 0; i < Vbound; i++){
      m = os[e].get_norb() + norb*n_bound[i][nrdim];
      res(m,m) += os[e].get_en();
    }
  }

  return res;
}
