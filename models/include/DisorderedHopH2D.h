#ifndef DISORDEREDHOPH2D_H
#define DISORDEREDHOPH2D_H

#include "TBDisorderedH2D.h"

using namespace std;

class DisorderedHopH2D : public TBDisorderedH2D{
  public:
    DisorderedHopH2D(TBModel model);
    ~DisorderedHopH2D();

    void generateDisorder();
    void setWeight(double w){this->w = w;}
    void setDisType(int type);

  private:

    int disType;
    void setIntraDisHop();
    double w;
};

#endif
