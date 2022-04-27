#ifndef TBDISORDEREDH_H
#define TBDISORDEREDH_H

#include "TBCleanH.h"

using namespace arma;
using namespace std;

class TBDisordered : public TBCleanH{
  cx_mat H(double * k);
};

#endif
