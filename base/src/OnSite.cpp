#include "OnSite.h"

using namespace std;

OnSite::OnSite(int nOrb, complex<double> en){
  this->nOrb = nOrb;
  this->en = en;
}

OnSite::OnSite(const OnSite & copy){
  nOrb = copy.nOrb;
  en = copy.en;
}

OnSite::~OnSite(){
}
