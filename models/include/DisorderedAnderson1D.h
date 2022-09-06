#ifndef DISORDEREDANDERSON1D_H
#define DISORDEREDANDERSON1D_H

#include "TBDisorderedH1D.h"

using namespace std;

class DisorderedAnderson1D : public TBDisorderedH1D{
  public:
    DisorderedAnderson1D(TBModel model);
    ~DisorderedAnderson1D();

    void generateDisorder();
    void setWeight(double w){this->w = w;}
    void setDisType(int type);

  private:

    int disType;
    void setIntraDisHop();
    double w;
};

#endif
