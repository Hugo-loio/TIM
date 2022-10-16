#ifndef DISORDEREDANDERSON3D_H
#define DISORDEREDANDERSON3D_H

#include "TBDisorderedH3D.h"

using namespace std;

class DisorderedAnderson3D : public TBDisorderedH3D{
  public:
    DisorderedAnderson3D(TBModel model);
    ~DisorderedAnderson3D();

    void generateDisorder();
    void setWeight(double w){this->w = w;}
    void setDisType(int type);

  private:

    int disType = 0;
    void setIntraDisHop();
    double w;
};

#endif
