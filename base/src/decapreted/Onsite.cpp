#include "Onsite.h"

using namespace std;

Onsite::Onsite(int norb, complex<double> en){
  this->norb = norb;
  this->en = en;
}

Onsite::Onsite(const Onsite & copy){
  norb = copy.norb;
  en = copy.en;
}

Onsite::~Onsite(){
}


