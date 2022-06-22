#include <iostream>
#include <armadillo>
#include "OData.h"
#include "SSH.h"
#include "SSH2D.h"
#include "BBH2D.h"
#include "BBH3D.h"
#include "SOTAI.h"

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
  bbh2D.setInterHop(2);
  bbh2D.getBands(argv[0], "EnergyBandsBBH2D_inter2.dat", 100, 100);

  bbh2D.setInterHop(1.5);
  bbh2D.getBands(argv[0], "EnergyBandsBBH2D_inter1.5.dat", 100, 100);

  bbh2D.setInterHop(1);
  bbh2D.getBands(argv[0], "EnergyBandsBBH2D_inter1.dat", 100, 100);

  bbh2D.setInterHop(0.5);
  bbh2D.getBands(argv[0], "EnergyBandsBBH2D_inter0.5.dat", 100, 100);

  cout << "3D BBH" << endl;
  BBH3D bbh3D(1,2);
  bbh3D.getBands(argv[0], "EnergyBandsBBH3D_inter2.dat");
  bbh3D.setInterHop(1);
  bbh3D.getBands(argv[0], "EnergyBandsBBH3D_inter1.dat");
  bbh3D.setInterHop(0.5);
  bbh3D.getBands(argv[0], "EnergyBandsBBH3D_inter0.5.dat");

  cout << "SOTAI" << endl;
  SOTAI sotai(2);
  sotai.getBands(argv[0], "EnergyBandsSOTAI_m2.dat", 100, 100);
  sotai.setM(1.01);
  sotai.getBands(argv[0], "EnergyBandsSOTAI_m1.01.dat", 100, 100);
  sotai.setM(1);
  sotai.getBands(argv[0], "EnergyBandsSOTAI_m1.dat", 100, 100);
  sotai.setM(0.99);
  sotai.getBands(argv[0], "EnergyBandsSOTAI_m0.99.dat", 100, 100);
  sotai.setM(0.5);
  sotai.getBands(argv[0], "EnergyBandsSOTAI_m0.5.dat", 100, 100);

  return 0;
}
