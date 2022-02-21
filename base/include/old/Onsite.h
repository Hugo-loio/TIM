#ifndef ONSITEOLD_H
#define ONSITEOLD_H

#include <complex>

using namespace std;

class Onsite{
  public:
    //On-site energy en at orbital norb
    Onsite(int norb, complex<double> en);
    Onsite(const Onsite & copy);
    ~Onsite();
    //Max distance between cells in a specific direction
    const int & get_norb(){return norb;};
    const complex<double> & get_en(){return en;};

  private:
    int norb;
    complex<double> en;
};

#endif
