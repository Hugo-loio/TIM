#ifndef DISORDEREDSSH_H
#define DISORDEREDSSH_H

#include "DisorderedHopH1D.h"

using namespace std;

class DisorderedSSH{
  public:
    DisorderedSSH(double t1 = 0.5, double t2 = 1);
    ~DisorderedSSH();

    void setIntraHop(double);
    void setInterHop(double);
    void setSize(int * l){ham->setSize(l);}
    void setW(double);

    void generateDisorder();

    double ipr(int nStates);
    double polarization();
    //double entanglement(int type);

    void test(char * argv0);

  private:
    TBModel * model;
    DisorderedHopH1D * ham;
};

#endif
