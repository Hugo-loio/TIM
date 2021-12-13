#ifndef SSH_1D_h
#define SSH_1D_h

#include "tbhop.h"
#include <armadillo>
using namespace arma;
using namespace std;

class SSH_1D : public tbhop::Kloop{


  public:
    SSH_1D(double t1 = 1, double t2 = 2); //t1 = intracell hopping , t2 = intercell hopping
    ~SSH_1D();

    void set_intrahop(double);
    void set_interhop(double);

    double get_interhop();
    double get_intrahop();

    cx_mat kH(double k); //Get hamiltonian in k space

    cx_mat hloop(double k){return kH(k);};

    vec eigenvalk(double k);
    cx_mat eigenveck(double k);

    double berry_kspace(int N); //N number of k intervals

    void set_rH(int L, bool closed); //Set hamiltonian in real space with L unit cells, closed = true(false) means periodic boundary conditions apply(don't apply). 
    sp_mat & get_rH(){return rH;};
    //double berry_twisted(int M, int N);  M is number of occupied states and N is number of k intervals

  private:
    double t1;
    double t2;

    sp_mat rH;
    int L;
    bool closed;
};

#endif
