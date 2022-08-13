#include "ParallelMPI.h"

ParallelMPI::ParallelMPI(int * argc, char *** argv){
  MPI_Init(argc, argv);
  MPI_Comm_size(MPI_COMM_WORLD, & nProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, & rank);
}

ParallelMPI::~ParallelMPI(){
  MPI_Finalize();
}
