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
    void setFile(char * argv0, string fileName); 
    vector<vector<double>> getResult(){return res;}

  private:
    bool printToFile;
    string fileName;
    int nThreads;
    vector<vector<double>> paramList;
    void (*job) (vector<double> &, vector<double>);
    OData * out;
    vector<vector<double>> res;
    chrono::high_resolution_clock::time_point tStart;

    void printTime();
};

#endif
