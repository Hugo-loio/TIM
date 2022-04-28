#ifndef HOPOLD_H
#define HOPOLD_H

#include <complex>

using namespace std;

class Hop{
  public:
    //Hopping hop between nOrb1 in home cell and nOrb2 in neighbouring cell defined by the nDim dimentional vector n 
    Hop(int nOrb1, int nOrb2, int * n, complex<double> hop, int nDim);
    Hop(const Hop & copy);
    ~Hop();

    //Max distance between cells in a specific direction
    int getMaxN(); 

    const int & getNOrb1(){return nOrb1;};
    const int & getNOrb2(){return nOrb2;};
    const int * getN(){return n;};
    const int & getN(int i){return n[i];};
    complex<double> getHop(){return hop;};
    void setHop(complex<double> hop){this->hop = hop;};

  private:
    int nOrb1;
    int nOrb2;
    int nDim;
    complex<double> hop;
    int * n;
};

#endif
