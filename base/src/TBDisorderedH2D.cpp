#include "TBDisorderedH2D.h"

TBDisorderedH2D::TBDisorderedH2D(TBModel model) : TBCleanH(model){
  nDisHop = 0;
  nDisOnSite = 0;
  isDisHop = new bool[nHop];
  isDisOnSite = new bool[nOnSite];
  indexDisHop = new int[nHop];
  indexDisOnSite = new int[nOnSite];
}

TBDisorderedH2D::~TBDisorderedH2D(){
  deleteDisArrays();
  delete[] isDisHop;
  delete[] isDisOnSite;
  delete[] indexDisHop;
  delete[] indexDisOnSite;
}

void TBDisorderedH2D::setDisHop(int nDisHop, int * hop){
  deleteDisArrays();
  for(int i = 0; i < nHop; i++){
    isDisHop[i] = false;
    indexDisHop[i] = 0;
  }
  for(int i = 0; i < nDisHop; i++){
    isDisHop[hop[i]] = true;
    indexDisHop[hop[i]] = i;
  }
  this->nDisHop = nDisHop;
}

void TBDisorderedH2D::setDisOnSite(int nDisOnSite, int * onSite){
  deleteDisArrays();
  for(int i = 0; i < nOnSite; i++){
    isDisOnSite[i] = false;
    indexDisOnSite[i] = 0;
  }
  for(int i = 0; i < nDisOnSite; i++){
    isDisOnSite[onSite[i]] = true;
    indexDisOnSite[onSite[i]] = i;
  }
  this->nDisOnSite = nDisOnSite;
}

void TBDisorderedH2D::setSize(int * l){
  deleteDisArrays();
  TBCleanH::setSize(l);
}

void TBDisorderedH2D::deleteDisArrays(){
  if(disHop != NULL && disOnSite != NULL){
    int e;
    for(int i = 0; i < nDisHop; i++){
      for(e = 0; e < l[0]; e++){
	delete[] disHop[i][e];
      }
      delete[] disHop[i];
    }
    for(int i = 0; i < nDisOnSite; i++){
      for(e = 0; e < l[0]; e++){
	delete[] disOnSite[i][e];
      }
      delete[] disOnSite[i];
    }
    delete[] disHop;
    delete[] disOnSite;
    disHop = NULL;
    disOnSite = NULL;
  }
}

void TBDisorderedH2D::createDisArrays(){
  disHop = new complex<double> ** [nDisHop];
  int e;
  for(int i = 0; i < nDisHop; i++){
    disHop[i] = new complex<double> * [l[0]];
    for(e = 0; e < l[1]; e++){
      disHop[i][e] = new complex<double> [l[1]];
    }
  }

  disOnSite = new complex<double> ** [nDisOnSite];
  for(int i = 0; i < nDisOnSite; i++){
    disOnSite[i] = new complex<double> * [l[0]];
    for(e = 0; e < l[1]; e++){
      disOnSite[i][e] = new complex<double> [l[1]];
    }
  }
}

cx_mat TBDisorderedH2D::H(double * k){
  if(!isUpdated){
    calcAux();
    isUpdated = true;
  }

  int size; 
  if(nRDim != 0){
    size = nOrb*lAccum[nRDim-1];
  }
  else{
    size = nOrb;
  }
  cx_mat res(size, size, fill::zeros);

  if(nRDim != 2 || nDim != 2){
    cout << "2D disordered class doesn't have a 2D TB model in real space. Expect things to go wrong, not my problem." << endl;
  }

  complex<double> t;
  int j,n,p,q,r;
  complex<double> phase;
  complex<double> ii(0,1);

  int * i = new int[nRDim + 1];
  int * b = new int[nRDim + 1];
  int * start = new int[nRDim];
  int * end = new int[nRDim + 1];
  end[nRDim] = 1;

  //Hopping terms
  for(int e = 0; e < nHop; e++){
    t = model.getHop(e).getHop();

    //Bulk lattice
    for(j = 0; j <= nRDim; j++){
      i[j] = startHopUCBulk[e][j];
    }

    n = flatten2(model.getHop(e).getNOrb1(), i);

    if(isDisHop[e] == false){
      while(i[nRDim] == 0){
	res(n, n + nHopBulk[e]) += t;

	i[0] += 1;
	n += incNHopBulk[e][0];
	p = 0;
	while(i[p] > endHopUCBulk[e][p]){
	  i[p] = startHopUCBulk[e][p];
	  i[++p] += 1;
	  n += incNHopBulk[e][p];
	}
      }
    }
    else{
      while(i[nRDim] == 0){
	res(n, n + nHopBulk[e]) += disHop[indexDisHop[e]][i[rOrder[0]]][i[rOrder[1]]];

	i[0] += 1;
	n += incNHopBulk[e][0];
	p = 0;
	while(i[p] > endHopUCBulk[e][p]){
	  i[p] = startHopUCBulk[e][p];
	  i[++p] += 1;
	  n += incNHopBulk[e][p];
	}
      }
    }

    //Boundary lattice
    for(j = 0; j < nRDim + 1; j++){
      b[j] = 0; //0 if not in boundary, 1 if in boundary for direction rIndex[j]
    }

    b[0]++;
    if(nRDim != 0){
      r = 0;
      while(b[r] > 1 || model.getHop(e).getN(rIndex[r]) == 0 || bC[rIndex[r]] == 0){
	b[r] = 0;
	b[++r]++;
	if(r == nRDim){
	  break;
	}
      }
    }

    while(b[nRDim] == 0){

      q = 0;
      phase = 1;
      for(j = 0; j < nRDim; j++){
	if(b[j] == 0){
	  start[j] = startHopUCBulk[e][j];
	  end[j] = endHopUCBulk[e][j];
	}
	else{
	  start[j] = startHopUCBound[e][j];
	  end[j] = endHopUCBound[e][j];
	  q += pow(2,rIndex[j]);
	  phase *= exp(-ii*theta[rIndex[j]]);
	}
	i[j] = start[j];
      }
      i[nRDim] = 0;

      t = model.getHop(e).getHop()*phase;

      n = flatten2(model.getHop(e).getNOrb1(), i);

      if(isDisHop[e] == false){
	while(i[nRDim] == 0){
	  res(n, n + nHopBound[e][q]) += t;

	  i[0] += 1;
	  n += incNHopBound[e][q][0];
	  p = 0;
	  while(i[p] > end[p]){
	    i[p] = start[p];
	    i[++p] += 1;
	    n += incNHopBound[e][q][p];
	  }
	}
      }
      else{
	while(i[nRDim] == 0){
	  res(n, n + nHopBound[e][q]) += disHop[indexDisHop[e]][i[rOrder[0]]][i[rOrder[1]]];

	  i[0] += 1;
	  n += incNHopBound[e][q][0];
	  p = 0;
	  while(i[p] > end[p]){
	    i[p] = start[p];
	    i[++p] += 1;
	    n += incNHopBound[e][q][p];
	  }
	}
      }


      b[0]++;
      r = 0;
      while(b[r] > 1 || model.getHop(e).getN(rIndex[r]) == 0 || bC[rIndex[r]] == 0){
	b[r] = 0;
	b[++r]++;
	if(r == nRDim){
	  break;
	}
      }
    }

  }
  res = res + res.t();

  //On-site terms

  for(int e = 0; e < nRDim; e++){
    end[e] = l[rIndex[e]] - 1;
  }

  for(int e = 0; e < nOnSite; e++){

    for(j = 0; j < nRDim; j++){
      i[j] = 0;
    }
    n = flatten2(model.getOnSite(e).getNOrb(), i);
    i[nRDim] = 0;

    if(isDisOnSite[e] == false){
      while(i[nRDim] == 0){
	res(n, n) += model.getOnSite(e).getEn();

	i[0] += 1;
	n += incNOnSite[e][0];
	p = 0;
	while(i[p] > end[p]){
	  i[p] = 0;
	  i[++p] += 1;
	  n += incNOnSite[e][p];
	}
      }
    }
    else{
      while(i[nRDim] == 0){
	res(n, n) += model.getOnSite(e).getEn();

	i[0] += 1;
	n += incNOnSite[e][0];
	p = 0;
	while(i[p] > end[p]){
	  i[p] = 0;
	  i[++p] += 1;
	  n += incNOnSite[e][p];
	}
      }
    }
  }

  delete[] i;
  delete[] b;
  delete[] start;
  delete[] end;

  return res;
}

sp_cx_mat TBDisorderedH2D::spH(double * k){
  if(!isUpdated){
    calcAux();
    isUpdated = true;
  }

  int size; 

  if(nRDim != 0){
    size = nOrb*lAccum[nRDim-1];
  }
  else{
    size = nOrb;
  }
  sp_cx_mat res(size, size);

  complex<double> t;
  int j,n,p,q,r;
  complex<double> phase;
  complex<double> ii(0,1);

  int * i = new int[nRDim + 1];
  int * b = new int[nRDim + 1];
  int * start = new int[nRDim];
  int * end = new int[nRDim + 1];
  end[nRDim] = 1;

  //Hopping terms
  for(int e = 0; e < nHop; e++){
    t = model.getHop(e).getHop();

    //Bulk lattice
    for(j = 0; j <= nRDim; j++){
      i[j] = startHopUCBulk[e][j];
    }

    n = flatten2(model.getHop(e).getNOrb1(), i);

    if(isDisHop[e] == false){
      while(i[nRDim] == 0){
	res(n, n + nHopBulk[e]) += t;

	i[0] += 1;
	n += incNHopBulk[e][0];
	p = 0;
	while(i[p] > endHopUCBulk[e][p]){
	  i[p] = startHopUCBulk[e][p];
	  i[++p] += 1;
	  n += incNHopBulk[e][p];
	}
      }
    }
    else{
      while(i[nRDim] == 0){
	res(n, n + nHopBulk[e]) += disHop[indexDisHop[e]][i[rOrder[0]]][i[rOrder[1]]];

	i[0] += 1;
	n += incNHopBulk[e][0];
	p = 0;
	while(i[p] > endHopUCBulk[e][p]){
	  i[p] = startHopUCBulk[e][p];
	  i[++p] += 1;
	  n += incNHopBulk[e][p];
	}
      }
    }

    //Boundary lattice
    for(j = 0; j < nRDim + 1; j++){
      b[j] = 0; //0 if not in boundary, 1 if in boundary for direction rIndex[j]
    }

    b[0]++;
    if(nRDim != 0){
      r = 0;
      while(b[r] > 1 || model.getHop(e).getN(rIndex[r]) == 0 || bC[rIndex[r]] == 0){
	b[r] = 0;
	b[++r]++;
	if(r == nRDim){
	  break;
	}
      }
    }

    while(b[nRDim] == 0){

      q = 0;
      phase = 1;
      for(j = 0; j < nRDim; j++){
	if(b[j] == 0){
	  start[j] = startHopUCBulk[e][j];
	  end[j] = endHopUCBulk[e][j];
	}
	else{
	  start[j] = startHopUCBound[e][j];
	  end[j] = endHopUCBound[e][j];
	  q += pow(2,rIndex[j]);
	  phase *= exp(-ii*theta[rIndex[j]]);
	}
	i[j] = start[j];
      }
      i[nRDim] = 0;

      t = model.getHop(e).getHop()*phase;

      n = flatten2(model.getHop(e).getNOrb1(), i);

      if(isDisHop[e] == false){
	while(i[nRDim] == 0){
	  res(n, n + nHopBound[e][q]) += t;

	  i[0] += 1;
	  n += incNHopBound[e][q][0];
	  p = 0;
	  while(i[p] > end[p]){
	    i[p] = start[p];
	    i[++p] += 1;
	    n += incNHopBound[e][q][p];
	  }
	}
      }
      else{
	while(i[nRDim] == 0){
	  res(n, n + nHopBound[e][q]) += disHop[indexDisHop[e]][i[rOrder[0]]][i[rOrder[1]]];

	  i[0] += 1;
	  n += incNHopBound[e][q][0];
	  p = 0;
	  while(i[p] > end[p]){
	    i[p] = start[p];
	    i[++p] += 1;
	    n += incNHopBound[e][q][p];
	  }
	}
      }


      b[0]++;
      r = 0;
      while(b[r] > 1 || model.getHop(e).getN(rIndex[r]) == 0 || bC[rIndex[r]] == 0){
	b[r] = 0;
	b[++r]++;
	if(r == nRDim){
	  break;
	}
      }
    }

  }
  res = res + res.t();

  //On-site terms

  for(int e = 0; e < nRDim; e++){
    end[e] = l[rIndex[e]] - 1;
  }

  for(int e = 0; e < nOnSite; e++){

    for(j = 0; j < nRDim; j++){
      i[j] = 0;
    }
    n = flatten2(model.getOnSite(e).getNOrb(), i);
    i[nRDim] = 0;

    if(isDisOnSite[e] == false){
      while(i[nRDim] == 0){
	res(n, n) += model.getOnSite(e).getEn();

	i[0] += 1;
	n += incNOnSite[e][0];
	p = 0;
	while(i[p] > end[p]){
	  i[p] = 0;
	  i[++p] += 1;
	  n += incNOnSite[e][p];
	}
      }
    }
    else{
      while(i[nRDim] == 0){
	res(n, n) += model.getOnSite(e).getEn();

	i[0] += 1;
	n += incNOnSite[e][0];
	p = 0;
	while(i[p] > end[p]){
	  i[p] = 0;
	  i[++p] += 1;
	  n += incNOnSite[e][p];
	}
      }
    }
  }

  delete[] i;
  delete[] b;
  delete[] start;
  delete[] end;

  return res;
}

cx_mat TBDisorderedH2D::blockH(int line, int col, double * k){
  if(!isUpdated){
    calcAux();
    isUpdated = true;
  }

  int size; 
  if(nRDim == 0){
    cout << "__PRETTY_FUNCTION__" << " can't divide systems in spatial layers, returning empty matrix" << endl;
    return cx_mat(0,0);
  }
  else{ 
    int nL = nuAccum[0]*nu[rIndex[0]];
    if(bDim == 0){
      size = nOrb/nL;
    }
    else{
      size = l[rIndex[0]]*(nOrb/nuAccum[0]);
      for(int i = 1; i < bDim; i++){
	size *= l[rIndex[i]]*nu[rIndex[i]];
      }
    }
    if(size*(col+1) > nOrb*lAccum[nRDim-1] || size*(line+1) > nOrb*lAccum[nRDim - 1]){
      return cx_mat(size,size, fill::zeros);
    }
  }
  cx_mat res(size, size, fill::zeros);

  int * startNL = new int[nRDim - bDim];
  int * startL = new int[nRDim - bDim];
  int * diffNL = new int[nRDim - bDim];
  int * diffN = new int[nRDim - bDim];
  int * diffL = new int[nRDim - bDim];

  bool transpose = false;
  if(col < line){
    int temp = col;
    col = line;
    line = temp;
    transpose = true;
  }

  diffNL[0] = (col-line) % (l[rIndex[bDim]]*nu[rIndex[bDim]]);
  startNL[0] = line % (l[rIndex[nDim]]*nu[rIndex[bDim]]);
  for(int i = 1; i < nRDim-bDim; i++){
    diffNL[i] = (col - line - diffNL[i-1])/((l[rIndex[bDim + i -1]]*nu[rIndex[bDim + i - 1]]));
    startNL[i] = line/((l[rIndex[bDim + i -1]]*nu[rIndex[bDim + i - 1]]));
    if(bDim + i < nRDim){
      diffNL[i] = diffNL[i] % (l[rIndex[bDim + i]]*nu[rIndex[bDim + i]]);
      startNL[i] = startNL[i] % (l[rIndex[bDim + i]]*nu[rIndex[bDim + i]]);
    }
  }

  for(int i = 0; i < nRDim - bDim; i++){
    startL[i] = startNL[i] % nu[rIndex[bDim + i]];
    diffN[i] = (startL[i] + diffNL[i])/nu[rIndex[bDim + i]];
    diffL[i] = diffNL[i] % nu[rIndex[bDim + i]];
    if((startL[i] + diffL[i]) % nu[rIndex[bDim + i]] < startL[i]){
      diffL[i] = -diffL[i];
    }
  }

  bool isHop;
  bool isHopHerm;
  vector<int> h;
  vector<int> hHerm; //Conjugate transpose the matrix comming from these
  int n, l1, l2;
  for(int e = 0; e < nHop; e++){
    isHop = true;
    isHopHerm = true;
    for(int i = 0; i < nRDim - bDim; i++){
      n = model.getHop(e).getN(rIndex[bDim + i]);
      l2 = lOrb[model.getHop(e).getNOrb2()][rIndex[bDim + i]];
      l1 = lOrb[model.getHop(e).getNOrb1()][rIndex[bDim + i]];
      if(n != diffN[i] || (l2 - l1) != diffL[i] || l1 != startL[i]){
	isHop = false;
      }
      if(n != - diffN[i] || (l2 - l1) != - diffL[i] || l2 != startL[i]){
	isHopHerm = false;
      }
    }
    if(isHop == true){
      h.push_back(e);
    }
    else if(isHopHerm == true){
      hHerm.push_back(e);
    }
  }

  int sub1 = (col-line)*size;
  int sub2 = startL[0]*size;
  for(int i = 1; i < nRDim - bDim; i++){
    sub2 += startL[i]*size*l[rIndex[bDim + i - 1]]*nu[rIndex[bDim + i - 1]];
  }

  vector<int> os;
  bool isOS;
  if(sub1 == 0){
    for(int e = 0; e < nOrb; e++){
      isOS = true;
      for(int i = 0; i < nRDim - bDim; i++){
	if(lOrb[e][rIndex[i + bDim]] != startL[i]){
	  isOS = false;
	}
      }
      if(isOS){
	os.push_back(e);
      }
    }
  }

  delete diffL;
  delete diffN;
  delete diffNL;
  delete startNL;
  delete startL;

  //TODO: until here some things can be precalculated to increase efficiency

  if(nRDim != 2 || nDim != 2){
    cout << "2D disordered class doesn't have a 2D TB model in real space. Expect things to go wrong, not my problem." << endl;
  }


  complex<double> t;
  int j,p,q,r;
  complex<double> phase;
  complex<double> ii(0,1);

  int * i = new int[nRDim + 1];
  int * b = new int[bDim + 1];
  int * start = new int[bDim];
  int * end = new int[bDim + 1];
  end[bDim] = 1;

  //Hopping terms (hermitian conjugate)
  for(int m = 0; m < hHerm.size(); m++){
    int e = hHerm[m];

    t = model.getHop(e).getHop();

    //Bulk lattice
    for(j = 0; j < bDim; j++){
      i[j] = startHopUCBulk[e][j];
      end[j] = endHopUCBulk[e][j];
    }
    for(j = bDim; j <= nRDim; j++){
      i[j] = 0;
    }

    n = flatten2(model.getHop(e).getNOrb1(), i) - sub2;
    i[bDim] = startHopUCBulk[e][nRDim];

    if(isDisHop[e] == false){
      while(i[bDim] == 0){
	res(n - sub1, n + nHopBulk[e]) += t;

	i[0] += 1;
	n += incNHopBulk[e][0];
	p = 0;
	while(i[p] > end[p]){
	  i[p] = startHopUCBulk[e][p];
	  i[++p] += 1;
	  n += incNHopBulk[e][p];
	}
      }
    }
    else{
      while(i[bDim] == 0){
	//TODO
	res(n - sub1, n + nHopBulk[e]) += disHop[indexDisHop[e]][i[rOrder[0]]][i[rOrder[1]]];

	i[0] += 1;
	n += incNHopBulk[e][0];
	p = 0;
	while(i[p] > end[p]){
	  i[p] = startHopUCBulk[e][p];
	  i[++p] += 1;
	  n += incNHopBulk[e][p];
	}
      }
    }

    //Boundary lattice
    for(j = 0; j < bDim + 1; j++){
      b[j] = 0; //0 if not in boundary, 1 if in boundary for direction rIndex[j]
    }

    b[0]++;
    if(bDim != 0){
      r = 0;
      while(b[r] > 1 || model.getHop(e).getN(rIndex[r]) == 0 || bC[rIndex[r]] == 0){
	b[r] = 0;
	b[++r]++;
	if(r == bDim){
	  break;
	}
      }
    }

    while(b[bDim] == 0){

      q = 0;
      phase = 1;
      for(j = 0; j < bDim; j++){
	if(b[j] == 0){
	  start[j] = startHopUCBulk[e][j];
	  end[j] = endHopUCBulk[e][j];
	}
	else{
	  start[j] = startHopUCBound[e][j];
	  end[j] = endHopUCBound[e][j];
	  q += pow(2,rIndex[j]);
	  phase *= exp(-ii*theta[rIndex[j]]);
	}
	i[j] = start[j];
      }
      for(j = bDim; j <= nRDim; j++){
	i[j] = 0;
      }

      t = model.getHop(e).getHop()*phase;

      n = flatten2(model.getHop(e).getNOrb1(), i) - sub2;

      while(i[bDim] == 0){

	res(n - sub1, n + nHopBound[e][q]) += t;

	i[0] += 1;
	n += incNHopBound[e][q][0];
	p = 0;
	while(i[p] > end[p]){
	  i[p] = start[p];
	  i[++p] += 1;
	  n += incNHopBound[e][q][p];
	}
      }


      b[0]++;
      r = 0;
      while(b[r] > 1 || model.getHop(e).getN(rIndex[r]) == 0 || bC[rIndex[r]] == 0){
	b[r] = 0;
	b[++r]++;
	if(r == bDim){
	  break;
	}
      }
    }

  }

  res = res.t();

  //Hopping terms 
  for(int m = 0; m < h.size(); m++){
    int e = h[m];

    t = model.getHop(e).getHop();

    //Bulk lattice
    for(j = 0; j < bDim; j++){
      i[j] = startHopUCBulk[e][j];
      end[j] = endHopUCBulk[e][j];
    }
    for(j = bDim; j <= nRDim; j++){
      i[j] = 0;
    }

    n = flatten2(model.getHop(e).getNOrb1(), i) - sub2;
    i[bDim] = startHopUCBulk[e][nRDim];

    while(i[bDim] == 0){
      res(n, n + nHopBulk[e] - sub1) += t;

      i[0] += 1;
      n += incNHopBulk[e][0];
      p = 0;
      while(i[p] > end[p]){
	i[p] = startHopUCBulk[e][p];
	i[++p] += 1;
	n += incNHopBulk[e][p];
	if(p == bDim){
	  break;
	}
      }
    }

    //Boundary lattice
    for(j = 0; j < bDim + 1; j++){
      b[j] = 0; //0 if not in boundary, 1 if in boundary for direction rIndex[j]
    }

    b[0]++;
    if(bDim != 0){
      r = 0;
      while(b[r] > 1 || model.getHop(e).getN(rIndex[r]) == 0 || bC[rIndex[r]] == 0){
	b[r] = 0;
	b[++r]++;
	if(r == bDim){
	  break;
	}
      }
    }

    while(b[bDim] == 0){

      q = 0;
      phase = 1;
      for(j = 0; j < bDim; j++){
	if(b[j] == 0){
	  start[j] = startHopUCBulk[e][j];
	  end[j] = endHopUCBulk[e][j];
	}
	else{
	  start[j] = startHopUCBound[e][j];
	  end[j] = endHopUCBound[e][j];
	  q += pow(2,rIndex[j]);
	  phase *= exp(-ii*theta[rIndex[j]]);
	}
	i[j] = start[j];
      }
      for(j = bDim; j <= nRDim; j++){
	i[j] = 0;
      }

      t = model.getHop(e).getHop()*phase;

      n = flatten2(model.getHop(e).getNOrb1(), i) - sub2;

      while(i[bDim] == 0){
	res(n, n + nHopBound[e][q] - sub1) += t;

	i[0] += 1;
	n += incNHopBound[e][q][0];
	p = 0;
	while(i[p] > end[p]){
	  i[p] = start[p];
	  i[++p] += 1;
	  n += incNHopBound[e][q][p];
	}
      }


      b[0]++;
      r = 0;
      while(b[r] > 1 || model.getHop(e).getN(rIndex[r]) == 0 || bC[rIndex[r]] == 0){
	b[r] = 0;
	b[++r]++;
	if(r == bDim){
	  break;
	}
      }
    }

  }

  if(transpose){
    res = res.t();
  }
  if(sub1 == 0){
    res = res + res.t();
  }

  //On-site terms

  for(int e = 0; e < bDim; e++){
    end[e] = l[rIndex[e]] - 1;
  }

  for(int m = 0; m < os.size(); m++){
    int e = os[m];

    for(j = 0; j < nRDim; j++){
      i[j] = 0;
    }
    n = flatten2(model.getOnSite(e).getNOrb(), i) - sub2;
    i[nRDim] = 0;

    while(i[bDim] == 0){
      res(n, n) += model.getOnSite(e).getEn();

      i[0] += 1;
      n += incNOnSite[e][0];
      p = 0;
      while(i[p] > end[p]){
	i[p] = 0;
	i[++p] += 1;
	n += incNOnSite[e][p];
      }
    }
  }

  delete[] i;
  delete[] b;
  delete[] start;
  delete[] end;

  return res;
}

