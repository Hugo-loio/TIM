#ifndef PARALLELMPI_H
#define PARALLELMPI_H

#include "mpi.h"

using namespace std;

class ParallelMPI{
  public:
    ParallelMPI(int * argc, char *** argv);
    ~ParallelMPI();
  private:
    int rank;
    int nProcs;
};
#endif
