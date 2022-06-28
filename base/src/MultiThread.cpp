#include "MultiThread.h"
#include <iostream>
#include <iomanip>
#include <thread>

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
  out = new OData(argv0, fileName);
  printToFile = true;
}

void MultiThread::run(){
  int nJobs = paramList.size();

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

  cout << "0 %\r";
  cout.flush();
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
	cout << setprecision(4) << 100*(double)count/(double)nJobs << " %\r";
	cout.flush();
      }
    }
  }
  cout << "100 %" << endl;
}
