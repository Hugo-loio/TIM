#ifndef DISORDEREDSSH_H
#define DISORDEREDSSH_H

#include "DisorderedHopH1D.h"
#include "Entanglement.h"

using namespace std;

class DisorderedSSH{
  public:
    DisorderedSSH(double t1 = 0.5, double t2 = 1, double delta = 0);
    ~DisorderedSSH();

    void setIntraHop(double);
    void setInterHop(double);
    void setSize(int * l);
    void setW(double);
    void setLayers(bool *);

    void generateDisorder();

    double ipr(int nStates, double en = 0);
    double polarization();
    double entanglement(int type);

    cx_mat getHam();
    void test(char * argv0);

  private:
    TBModel * model;
    DisorderedHopH1D * ham;
    Entanglement * ent = NULL; 
    bool diagEnt = true;
};

#endif
