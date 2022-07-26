#ifndef DISORDEREDHOPH1D_H
#define DISORDEREDHOPH1D_H

#include "TBDisorderedH1D.h"

using namespace std;

class DisorderedHopH1D : public TBDisorderedH1D{
  public:
    DisorderedHopH1D(TBModel model);
    ~DisorderedHopH1D();

    void generateDisorder();
    void setWeight(double w){this->w = w;}
    void setDisType(int type);

  private:

    int disType;
    void setIntraDisHop();
    double w;
};

#endif
