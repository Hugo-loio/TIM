#ifndef HOP_H
#define HOP_H

#include <complex>

using namespace std;

class Hop{
  public:
    //Hopping hop between norb1 in home cell and orb2 in neighbouring cell defined by vector the ndim dimentional n 
    Hop(int norb1, int norb2, int * n, complex<double> hop, int ndim);
    Hop(const Hop & copy);
    ~Hop();
    //Max distance between cells in a specific direction
    int get_maxn(); 

    const int & get_norb1(){return norb1;};
    const int & get_norb2(){return norb2;};
    const int * get_n(){return n;};
    const complex<double> & get_hop(){return hop;};

  private:
    int norb1;
    int norb2;
    int ndim;
    complex<double> hop;
    int * n;
};

#endif
