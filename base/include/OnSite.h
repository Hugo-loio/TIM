#ifndef ONSITE_H
#define ONSITE_H

#include <complex>

using namespace std;

class OnSite{
  public:
    //On-site energy en at orbital nOrb
    OnSite(int nOrb, complex<double> en);
    OnSite(const OnSite & copy);
    ~OnSite();

    const int & getNOrb(){return nOrb;};
    const complex<double> & getEn(){return en;};
    void setEn(complex<double> en){this->en = en;};

  private:
    int nOrb;
    complex<double> en;
};

#endif
