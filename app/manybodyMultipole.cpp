#include <iostream>
#include "SSH.h"

int main (int arc, char ** argv) {
  SSH ssh(2,1);

  cout << ssh.polarization(10) << endl;
  cout << ssh.berryPhase(10) << endl;

  ssh.setIntraHop(0.5);
  cout << ssh.polarization(10) << endl;
  cout << ssh.berryPhase(10) << endl;
}
