// basic file operations
#include <iostream>
#include <armadillo>
#include "OData.h"
#include "SSH.h"
#include "SSH2D.h"
#include "BBH2D.h"
#include "BBH3D.h"

using namespace std;
using namespace arma;

int main (int arc, char ** argv) {
  ofstream f;
  f.open("test.txt");

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
  bbh2D.getWannierBands(argv[0], "WannierBandsBBH2D_inter2_y.dat", 1); 
  bbh2D.getWannierBands(argv[0], "WannierBandsBBH2D_inter2_x.dat", 0); 

  bbh2D.setInterHop(1.5);
  bbh2D.getBands(argv[0], "EnergyBandsBBH2D_inter1.5.dat", 100, 100);
  bbh2D.getWannierBands(argv[0], "WannierBandsBBH2D_inter1.5_y.dat", 1); 
  bbh2D.getWannierBands(argv[0], "WannierBandsBBH2D_inter1.5_x.dat", 0); 

  bbh2D.setInterHop(3);
  bbh2D.getWannierBands(argv[0], "WannierBandsBBH2D_inter3_y.dat", 1); 
  bbh2D.getWannierBands(argv[0], "WannierBandsBBH2D_inter3_x.dat", 0); 

  bbh2D.setInterHop(0.9);
  bbh2D.getWannierBands(argv[0], "WannierBandsBBH2D_inter0.9_x.dat", 0); 
  bbh2D.getWannierBands(argv[0], "WannierBandsBBH2D_inter0.9_y.dat", 1); 

  bbh2D.setInterHop(0.99);
  bbh2D.getWannierBands(argv[0], "WannierBandsBBH2D_inter0.99_x.dat", 0); 
  bbh2D.getWannierBands(argv[0], "WannierBandsBBH2D_inter0.99_y.dat", 1); 

  bbh2D.setInterHop(1.01);
  bbh2D.getWannierBands(argv[0], "WannierBandsBBH2D_inter1.01_x.dat", 0); 
  bbh2D.getWannierBands(argv[0], "WannierBandsBBH2D_inter1.01_y.dat", 1); 

  bbh2D.setInterHop(1.1);
  bbh2D.getWannierBands(argv[0], "WannierBandsBBH2D_inter1.1_x.dat", 0); 
  bbh2D.getWannierBands(argv[0], "WannierBandsBBH2D_inter1.1_y.dat", 1); 

  bbh2D.setInterHop(1);
  bbh2D.getWannierBands(argv[0], "WannierBandsBBH2D_inter1_x.dat", 0); 
  bbh2D.getWannierBands(argv[0], "WannierBandsBBH2D_inter1_y.dat", 1); 
  bbh2D.getBands(argv[0], "EnergyBandsBBH2D_inter1.dat", 100, 100);

  bbh2D.setInterHop(0.5);
  bbh2D.getBands(argv[0], "EnergyBandsBBH2D_inter0.5.dat", 100, 100);
  bbh2D.getWannierBands(argv[0], "WannierBandsBBH2D_inter0.5_y.dat", 1); 
  bbh2D.getWannierBands(argv[0], "WannierBandsBBH2D_inter0.5_x.dat", 0); 

  cout << "3D BBH" << endl;
  BBH3D bbh3D(1,2);
  double k[3] = {M_PI/2, M_PI/2, M_PI/2};
  bbh3D.getBands(argv[0], "EnergyBandsBBH3D_inter2.dat");
  bbh3D.getWannierBands(argv[0], "WannierBandsBBH3D_inter2_x.dat", 0); 
  bbh3D.getWannierBands(argv[0], "WannierBandsBBH3D_inter2_y.dat", 1); 
  bbh3D.getWannierBands(argv[0], "WannierBandsBBH3D_inter2_z.dat", 2); 
  //f << bbh3D.getH(k) << endl;
  bbh3D.setInterHop(1);
  bbh3D.getBands(argv[0], "EnergyBandsBBH3D_inter1.dat");
  //bbh3D.getWannierBands(argv[0], "WannierBandsBBH3D_inter1_x.dat", 0); 
  //f << bbh3D.getH(k) << endl;
  bbh3D.setInterHop(0.5);
  bbh3D.getBands(argv[0], "EnergyBandsBBH3D_inter0.5.dat");
  //bbh3D.getWannierBands(argv[0], "WannierBandsBBH3D_inter0.5_x.dat", 0); 

  f.close();
  return 0;
}
