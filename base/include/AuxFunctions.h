#include "MultiThread.h"
#include <chrono>

void run(void (*job) (vector<double> &, vector<double>), vector<vector<double>> & paramList, int threadNumber, char * argv0, string fileName, int nSamples = 1, int version = 0, int part = 0);

void runSingleThread(void (*job) (vector<double> &, vector<double>), vector<vector<double>> & paramList, char * argv0, string fileName, int nSamples = 1, int version = 0, int part = 0);

string rmTrailZeros(string str);

double cot(double x);

void printElapsedTime(chrono::high_resolution_clock::time_point tStart);
