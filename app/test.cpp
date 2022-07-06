#include <iostream>
#include <time.h>
#include <armadillo>
#include "OData.h"
#include "DisorderedSOTAI.h"
#include "BBH2D.h"
#include "BBH3D.h"
#include <thread>
#include <cstdlib>

using namespace std;
using namespace arma;

void dummy_func(){
  sleep(1);
  cout << "Hello world" << endl;
}

int main (int arc, char ** argv) {
  DisorderedSOTAI sotai(1.1);
  sotai.test(argv[0]);
  /*
  int l[2] = {10,10};
  sotai.setSize(l);
  sotai.setW(3);
  for(int i = 0; i < 10; i++){
    sotai.generateDisorder();
    cout << sotai.getQuadrupoleManyBody() << endl;
  }
  */

  //unsigned int n = thread::hardware_concurrency();
  //cout << n << " concurrent threads are supported.\n";

  /*
     BBH3D bbh3(0.5,1,0.5);
     bbh3.test(argv[0]);
     BBH2D bbh2(0.5,1,0.5);
     bbh2.test(argv[0]);
     */

  return 0;
}
