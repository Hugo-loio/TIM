#include <iostream>
#include <fstream>
#include <armadillo>
#include "OData.h"
#include "SSH.h"
#include "SSH2D.h"
#include "BBH2D.h"

using namespace std;
using namespace arma;

int main (int arc, char ** argv) {
  ofstream f;
  f.open("test.txt");

  cout << "1D SSH" << endl;
  SSH ssh(1,2);
  cout << "SSH Berry phase inter 2: " << ssh.berryPhase(10) << endl;
  cout << "Supercell: " << ssh.berryPhaseSupercell(10,10) << endl;
  ssh.setInterHop(1);
  cout << "SSH Berry phase inter 1: " << ssh.berryPhase(10) << endl;
  cout << "Supercell: " << ssh.berryPhaseSupercell(10,10) << endl;
  ssh.setInterHop(0.5);
  cout << "SSH Berry phase inter 0.5: " << ssh.berryPhase(10) << endl;
  cout << "Supercell: " << ssh.berryPhaseSupercell(10,10) << endl;

  double k0[2] = {0, -M_PI};
  double k[2] ={M_PI/2, M_PI/2};
  int bC[2] = {2,1};
  int l[2] = {10,10};
  int l2[2] = {2,2};
  int l3[2] = {3,3};
  cout << "\n2D SSH" << endl;
  SSH2D ssh2D(1,2);
  cout << "2D SSH Berry phase inter 2: " << ssh2D.berryPhase(100,0,k0) << endl;
  cout << "Supercell in x, reciprocal y: " << ssh2D.berryPhaseSupercell(10,0,bC,l,k0) << endl;
  ssh2D.setInterHop(1.5);
  cout << "2D SSH Berry phase inter 1.5: " << ssh2D.berryPhase(100,0,k0) << endl;
  cout << "Supercell in x, reciprocal y: " << ssh2D.berryPhaseSupercell(10,0,bC,l,k0) << endl;
  ssh2D.setInterHop(1);
  cout << "2D SSH Berry phase inter 1: " << ssh2D.berryPhase(100,0,k0) << endl;
  ssh2D.setInterHop(0.5);
  cout << "2D SSH Berry phase inter 0.5: " << ssh2D.berryPhase(100,0,k0) << endl;
  cout << "Supercell in x, reciprocal y: " << ssh2D.berryPhaseSupercell(10,0,bC,l,k0) << endl;
  //cout << "Supercell: " << ssh2D.berryPhaseSupercell(10,0,10,10) << endl;

  cout << "\n2D BBH" << endl;
  BBH2D bbh2D(1,2);
  cout << "2D BBH Berry phase inter 2: " << bbh2D.berryPhase(100,0,k0) << endl;
  bbh2D.setInterHop(1);
  cout << "2D BBH Berry phase inter 1: " << bbh2D.berryPhase(100,0,k0) << endl;
  bbh2D.setInterHop(0.5);
  cout << "2D BBH Berry phase inter 0.5: " << bbh2D.berryPhase(100,0,k0) << endl;

  f.close();
  return 0;
}
