#ifndef HOP_H
#define HOP_H

class hop{
  public:
    //Hopping hop between norb1 in home cell and orb2 in neighbouring cell defined by vector the ndim dimentional n 
    hop(int norb1, int norb2, int * n, double hop, int ndim);
    ~hop();
    //Max distance between cells in a specific direction
    int get_maxn(); 
    //Getting hopping amplitude
    double get_hop(){return hop;};

  private:
    int norb1;
    int norb2;
    int ndim;
    int * n;
    double hop;
};

#endif
