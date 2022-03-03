// basic file operations
#include <iostream>
#include <armadillo>
#include "OData.h"
#include "SSH.h"

using namespace std;
using namespace arma;

int main (int arc, char ** argv) {
  SSH ssh(1,2);
  ssh.getBands(argv[0], "EnergyBandsSSH.dat", 100);
  return 0;
}
