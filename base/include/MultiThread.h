#ifndef MULTITHREAD_H
#define MULTITHREAD_H

#include <vector>
#include "OData.h"
#include <string>

using namespace std;

class MultiThread{
  public:
    MultiThread(void (*job) (vector<double> & res, vector<double> params), vector<vector<double>> paramList, int nThreads);
    ~MultiThread();
    void run();
    //void runSingleThread();

    void setFile(char * argv0, string fileName); 
    void setAuxFile(char * argv0, string fileName);
    //Samples per thread
    void setSamples(int nSamples);
    vector<vector<double>> getResult(){return res;}

  private:
    bool printToFile;
    bool existsAuxFile;
    string fileName;
    int nThreads;
    vector<vector<double>> paramList;
    void (*job) (vector<double> &, vector<double>);
    OData * out;
    vector<vector<double>> res;
    chrono::high_resolution_clock::time_point tStart;
    int nSamples;

    void printTime();
    //void update
};

#endif
