// basic file operations
#include <iostream>
#include <armadillo>
#include "OData.h"
#include "SSH.h"
#include "SSH2D.h"
#include "BBH2D.h"

using namespace std;
using namespace arma;

int main (int arc, char ** argv) {

  cout << "1D SSH" << endl;
  SSH ssh(1,2);
  cout << "SSH Berry phase inter 2: " << ssh.berryPhase(10) << endl;
  ssh.getBands(argv[0], "EnergyBandsSSH_inter2.dat", 100);
  ssh.setInterHop(1);
  cout << "SSH Berry phase inter 1: " << ssh.berryPhase(10) << endl;
  ssh.getBands(argv[0], "EnergyBandsSSH_inter1.dat", 100);
  ssh.setInterHop(0.5);
  cout << "SSH Berry phase inter 0.5: " << ssh.berryPhase(10) << endl;
  ssh.getBands(argv[0], "EnergyBandsSSH_inter0.5.dat", 100);

  double k0[2] = {0, 0};
  cout << "2D SSH" << endl;
  SSH2D ssh2D(1,2);
  cout << "2D SSH Berry phase inter 2: " << ssh2D.berryPhase(100,0,k0) << endl;
  ssh2D.getBands(argv[0], "EnergyBandsSSH2D_inter2.dat", 20, 20);
  ssh2D.setInterHop(1);
  cout << "2D SSH Berry phase inter 1: " << ssh2D.berryPhase(100,0,k0) << endl;
  ssh2D.getBands(argv[0], "EnergyBandsSSH2D_inter1.dat", 20, 20);
  ssh2D.setInterHop(0.5);
  cout << "2D SSH Berry phase inter 0.5: " << ssh2D.berryPhase(100,0,k0) << endl;
  ssh2D.getBands(argv[0], "EnergyBandsSSH2D_inter0.5.dat", 20, 20);

  cout << "2D BBH" << endl;
  BBH2D bbh2D(1,2);
  cout << "2D BBH Berry phase inter 2: " << bbh2D.berryPhase(100,0,k0) << endl;
  double avgBerry = 0;
  for(int i = 0; i < 100; i++){
    avgBerry += bbh2D.berryPhase(100,0,k0);
  }
  avgBerry /= 100;
  cout << "Average Berry Phase: " << avgBerry << endl;
  bbh2D.getBands(argv[0], "EnergyBandsBBH2D_inter2.dat", 20, 20);
  bbh2D.setInterHop(1);
  cout << "2D BBH Berry phase inter 1: " << bbh2D.berryPhase(100,0,k0) << endl;
  avgBerry = 0;
  for(int i = 0; i < 100; i++){
    avgBerry += bbh2D.berryPhase(100,0,k0);
  }
  avgBerry /= 100;
  cout << "Average Berry Phase: " << avgBerry << endl;
  bbh2D.getBands(argv[0], "EnergyBandsBBH2D_inter1.dat", 20, 20);
  bbh2D.setInterHop(0.5);
  cout << "2D BBH Berry phase inter 0.5: " << bbh2D.berryPhase(100,0,k0) << endl;
  avgBerry = 0;
  for(int i = 0; i < 100; i++){
    avgBerry += bbh2D.berryPhase(100,0,k0);
  }
  avgBerry /= 100;
  cout << "Average Berry Phase: " << avgBerry << endl;
  bbh2D.getBands(argv[0], "EnergyBandsBBH2D_inter0.5.dat", 20, 20);

  return 0;
}
