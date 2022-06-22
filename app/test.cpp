#include <iostream>
#include <armadillo>
#include "OData.h"
#include "DisorderFunctions.h"
#include "DisorderedSOTAI.h"
#include "BBH2D.h"
#include "BBH3D.h"
#include <thread>
#include <cstdlib>

using namespace std;
using namespace arma;

class A{
  public:
    void func(){
      overfunc();
    };
  protected:
    void overfunc(){
      cout << "class A" << endl;
    }
};

class B: public A{
  public:
    void func(){
      A::func();
    };
  protected:
    void overfunc(){
      cout << "class B" << endl;
    };
};

int main (int arc, char ** argv) {
  /*
     for(int i = 0; i < 100; i++){
     cout << uniform(0,1) << endl;
     }
     */

  DisorderedSOTAI sotai(1.1);
  sotai.setW(1);


  //OData o(argv[0], "testH.dat");
  //o.matrixWeights(h);

  //int l[2] = {2,10};
  //vec eigVal;
  //cx_mat eigVec;
  //eig_sym(eigVal, eigVec, sotai.getHam(l));
  //sotai.getQuadrupoleManyBody(l);
  //cout << sotai.getHam(l) << endl;
  //sotai.getTopInv(l);

  //unsigned int n = thread::hardware_concurrency();
  //cout << n << " concurrent threads are supported.\n";

  BBH3D bbh3(0.5,1,3);
  bbh3.test(argv[0]);

  BBH2D bbh2(0.7,1);
  int l[2] = {100,100};
  //cout << bbh2.getBoundPolarization(l,0) << endl;
  //cout << bbh2.getBoundPolarization(l,1) << endl;

  /*
  A a;
  B b;
  b.func();
  */
  return 0;
}
