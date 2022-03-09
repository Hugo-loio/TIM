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
  ssh.getBands(argv[0], "EnergyBandsSSH_inter2.dat", 100);
  ssh.setInterHop(1);
  ssh.getBands(argv[0], "EnergyBandsSSH_inter1.dat", 100);
  ssh.setInterHop(0.5);
  ssh.getBands(argv[0], "EnergyBandsSSH_inter0.5.dat", 100);

  cout << "2D SSH" << endl;
  SSH2D ssh2D(1,2);
  ssh2D.getBands(argv[0], "EnergyBandsSSH2D_inter2.dat", 100, 100);
  ssh2D.setInterHop(1);
  ssh2D.getBands(argv[0], "EnergyBandsSSH2D_inter1.dat", 100, 100);
  ssh2D.setInterHop(0.5);
  ssh2D.getBands(argv[0], "EnergyBandsSSH2D_inter0.5.dat", 100, 100);

  cout << "2D BBH" << endl;
  BBH2D bbh2D(1,2);
  bbh2D.getBands(argv[0], "EnergyBandsBBH2D_inter2.dat", 100, 100);
  bbh2D.getWannierBands(argv[0], "WannierBandsBBH2D_inter2_y.dat", 1, 10); 
  bbh2D.getWannierBands(argv[0], "WannierBandsBBH2D_inter2_x.dat", 0, 10); 

  bbh2D.setInterHop(1);
  bbh2D.getBands(argv[0], "EnergyBandsBBH2D_inter1.dat", 100, 100);

  bbh2D.setInterHop(0.5);
  bbh2D.getBands(argv[0], "EnergyBandsBBH2D_inter0.5.dat", 100, 100);
  bbh2D.getWannierBands(argv[0], "WannierBandsBBH2D_inter0.5_y.dat", 1, 10); 
  bbh2D.getWannierBands(argv[0], "WannierBandsBBH2D_inter0.5_x.dat", 0, 10); 

  return 0;
}
