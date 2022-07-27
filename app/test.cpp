#include <iostream>
#include <time.h>
#include <armadillo>
#include "OData.h"
#include "DisorderedSOTAI.h"
#include "DisorderedSSH.h"
#include "BBH2D.h"
#include "BBH3D.h"
#include <thread>
#include <cstdlib>

using namespace std;
using namespace arma;

int main (int arc, char ** argv) {
  /*
     DisorderedSOTAI sotai(1.1);
     int l[2] = {20,20};
     sotai.setSize(l);
     sotai.setW(0);
     for(int i = 0; i < 10; i++){
     sotai.generateDisorder();
     sotai.getIPR(10);
     }
  //sotai.getTMM(100, 0, 10);
  */

  DisorderedSSH ssh(0.5);
  int l[1] = {10};
  ssh.setSize(l);
  ssh.setW(0.5);
  ssh.generateDisorder();
  vec eigval;
  eig_sym(eigval, ssh.getHam());
  cout << eigval << endl;
  cout << ssh.getHam() << endl;

  return 0;
}
