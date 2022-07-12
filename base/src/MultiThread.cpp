#include "MultiThread.h"
#include <iostream>
#include <iomanip>
#include <thread>
#include <chrono>

MultiThread::MultiThread(void (*job) (vector<double> &, vector<double>), vector<vector<double>> paramList, int nThreads){
  printToFile = false;
  this->nThreads = nThreads;
  this->paramList = paramList;
  this->job = job;
}

MultiThread::~MultiThread(){
  delete out;
}

void MultiThread::setFile(char * argv0, string fileName){
  this->fileName = fileName;
  out = new OData(argv0, fileName);
  printToFile = true;
}

void MultiThread::run(){
  int nJobs = paramList.size();
  if(printToFile){
    cout << "Printing to file: " << fileName << endl;
  }
  cout << "Running " << nJobs << " jobs in " << nThreads << "threads." << endl;

  res.clear();
  for(int i = 0; i < nJobs; i++){
    res.push_back(vector<double>());
  }

  vector<thread> t; 
  if(nJobs < nThreads){
    nThreads = nJobs;
  }
  int * resIndex = new int[nThreads];
  for(int i = 0; i < nThreads; i++){
    t.push_back(thread(job, ref(res[i]), paramList[i]));
    resIndex[i] = i;
  }

  tStart = chrono::high_resolution_clock::now();
  cout << "0 \% at ";
  printTime();
  int count = 0;
  while(count < nJobs){
    for(int i = 0; i < nThreads; i++){
      if(t[i].joinable()){
	t[i].join();
	if(printToFile){
	  out->line(res[resIndex[i]]);
	}
	if(count + nThreads < nJobs){
	  t[i] = thread(job, ref(res[count + nThreads]), paramList[count + nThreads]);
	  resIndex[i] = count + nThreads;
	}
	count++;
	cout << setprecision(4) << 100*(double)count/(double)nJobs << " \% at thread " << i << " at ";
	printTime();
      }
    }
  }
  cout << "100 %" << endl;
  delete[] resIndex;
}

void MultiThread::printTime(){
  auto tNow = chrono::high_resolution_clock::now();
  auto dH = chrono::duration_cast<chrono::hours>(tNow-tStart);
  auto dM = chrono::duration_cast<chrono::minutes>(tNow-tStart);
  auto dS = chrono::duration_cast<chrono::seconds>(tNow-tStart);
  int h = dH.count();
  int m = dM.count() % 60;
  int s = dS.count() % 60;
  cout << h << ":" << m << ":" << s << endl;
}

