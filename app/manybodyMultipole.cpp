#include <iostream>
#include <iomanip>
#include "SSH.h"
#include "BBH2D.h"
#include "BBH3D.h"

int main (int arc, char ** argv) {
  //SSH
  SSH ssh(2,1);
  cout << "SSH" << endl;
  cout << "Trivial" << endl;
  cout << "Many-body: " << ssh.polarization(10) << " Berry phase: " <<  ssh.berryPhase(10) << endl;

  cout << "Topological" << endl;
  ssh.setIntraHop(0.5);
  cout << "Many-body: " << ssh.polarization(10) << " Berry phase: " <<  ssh.berryPhase(10) << endl;

  //BBH2D
  BBH2D bbh2D(2,1);
  int l2D[2] = {10,10};
  cout << "\nBBH2D" << endl;
  cout << "Trivial" << endl;
  cout << "Many-body: " << bbh2D.getQuadrupoleManyBody(l2D) << " Nested Wilson: " <<  bbh2D.getQuadrupoleNested(10,10) << endl;

  cout << "Topological" << endl;
  bbh2D.setIntraHop(0.5);
  cout << "Many-body: " << bbh2D.getQuadrupoleManyBody(l2D) << " Nested Wilson: " <<  bbh2D.getQuadrupoleNested(10,10) << endl;

  //BBH3D
  BBH3D bbh3D(2,1);
  int l3D[3] = {4,4,4};
  cout << "\nBBH2D" << endl;
  cout << "Trivial" << endl;
  cout << "Many-body: " << bbh3D.getOctupoleManyBody(l3D) << " Nested Wilson: " <<  bbh3D.getOctupoleNested(4,4,4) << endl;

  cout << "Topological" << endl;
  bbh3D.setIntraHop(0.5);
  cout << "Many-body: " << bbh3D.getOctupoleManyBody(l3D) << " Nested Wilson: " <<  bbh3D.getOctupoleNested(4,4,4) << endl;

}
