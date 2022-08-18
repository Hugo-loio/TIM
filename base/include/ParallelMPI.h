#ifndef PARALLELMPI_H
#define PARALLELMPI_H

#include "mpi.h"
#include "OData.h"

using namespace std;

class ParallelMPI{
  public:
    ParallelMPI(int * argc, char *** argv);
    ~ParallelMPI();

    void setFile(char * argv0, string fileName); 
    void setSamples(int nSamples);
    void setParamList(double ** paramList, int paramSize, int nParams);
    void setParamList(vector<vector<double>> paramList);
    void setJob(void (*job) (double * res, double * params), int resSize);
    void setPrintEachSamp(bool printEachSamp){this->printEachSamp = printEachSamp;};

    void run();

    vector<vector<double>> getResult(){return fullRes;}

  private:
    //General attributes
    int nSamples = 1;
    string fileName;
    OData * out;
    double ** paramList;
    int paramSize;
    int nParams;
    int resSize;
    void (*job) (double *, double *);
    chrono::high_resolution_clock::time_point tStart;
    vector<vector<double>> fullRes;
    bool printEachSamp = false;

    //MPI related attributes
    int rank;
    int nProcs;
    int jobCount = 0;

    //Checks 
    bool isJobSet = false;
    bool areParamsSet = false;
    bool isFileSet = false;

    void createType();
    void printTime();
    void printTimePerJob();
};
#endif
