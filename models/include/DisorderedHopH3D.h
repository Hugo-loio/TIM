#ifndef DISORDEREDHOPH3D_H
#define DISORDEREDHOPH3D_H

#include "TBDisorderedH3D.h"

using namespace std;

class DisorderedHopH3D : public TBDisorderedH3D{
  public:
    DisorderedHopH3D(TBModel model);
    ~DisorderedHopH3D();

    void generateDisorder();
    void setWeight(double w){this->w = w;}
    void setDisType(int type);

  private:

    int disType;
    void setIntraDisHop();
    double w;
};

#endif
