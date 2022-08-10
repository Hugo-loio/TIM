#include "MultiThread.h"
#include <iostream>
#include <iomanip>
#include <thread>
#include <chrono>

MultiThread::MultiThread(void (*job) (vector<double> &, vector<double>), vector<vector<double>> paramList, int nThreads){
  printToFile = false;
  existsAuxFile = false;
  this->nThreads = nThreads;
  this->paramList = paramList;
  this->job = job;
  nSamples = 1;
}

MultiThread::~MultiThread(){
  delete out;
}

void MultiThread::setFile(char * argv0, string fileName){
  this->fileName = fileName;
  out = new OData(argv0, fileName);
  printToFile = true;
}

void MultiThread::setSamples(int nSamples){
  this->nSamples = nSamples;
}

void MultiThread::run(){
  int nJobs = paramList.size();
  if(printToFile){
    cout << "Printing to file: " << fileName << endl;
  }
  if(existsAuxFile){
    cout << "Auxiliary file: " << auxFileName << endl;
  }
  cout << "Running " << nJobs << " jobs in " << nThreads << "threads." << endl;

  res.clear();
  for(int i = 0; i < nJobs*nSamples; i++){
    res.push_back(vector<double>());
  }

  vector<thread> t; 
  if(nJobs < nThreads){
    nThreads = nJobs;
  }
  int * resIndex = new int[nThreads];
  int * completSamp = new int[nJobs];
  int * startedSamp = new int[nJobs];
  for(int i = 0; i < nJobs; i++){
    completSamp[i] = 0;
    startedSamp[i] = 0;
  }
  int e;
  int countStart = 0;
  for(int i = 0; i < nThreads; i++){
    e = i/nSamples;
    t.push_back(thread(job, ref(res[i]), paramList[e]));
    resIndex[i] = i;
    startedSamp[e]++;
    if(startedSamp[e] == nSamples){
      countStart++;
    }
  }

  tStart = chrono::high_resolution_clock::now();
  cout << "0 \% at ";
  printTime();
  int countEnd = 0;
  int jobIndex = 0;
  vector<double> data;
  while(countEnd < nJobs){
    for(int i = 0; i < nThreads; i++){
      if(t[i].joinable()){
	t[i].join();
	jobIndex = resIndex[i]/nSamples;
	completSamp[jobIndex]++;
	//cout << countStart << " " << countEnd << " " << resIndex[i] << " " << jobIndex << " " <<  completSamp[jobIndex] <<  " " << nJobs << endl;
	if(completSamp[jobIndex] == nSamples){
	  if(printToFile){
	    data = paramList[jobIndex];
	    for(e = 0; e < nSamples; e++){
	      data.insert(data.end(), res[jobIndex*nSamples + e].begin(), res[jobIndex*nSamples + e].end());
	    }
	    out->line(data);
	  }
	  countEnd++;
	  cout << setprecision(4) << 100*(double)countEnd/(double)nJobs << " \% at thread " << i << " at ";
	  printTime();
	}
	if(countStart < nJobs){
	  e = (countStart)*nSamples + startedSamp[countStart];
	  t[i] = thread(job, ref(res[e]), paramList[countStart]);
	  resIndex[i] = e;
	  startedSamp[countStart]++;
	}
	if(startedSamp[countStart] == nSamples){
	  countStart++;
	}
      }
    }
  }
  cout << "100 %" << endl;
  delete[] resIndex;
  delete[] completSamp;
  delete[] startedSamp;
}

/*
   void MultiThread::runSingleThread(){
   int nJobs = paramList.size();
   if(printToFile){
   cout << "Printing to file: " << fileName << endl;
   }
   cout << "Running " << nJobs << " jobs in " << nThreads << "threads." << endl;

   res.clear();
   for(int i = 0; i < nJobs; i++){
   res.push_back(vector<double>());
   }

   tStart = chrono::high_resolution_clock::now();

   cout << "0 \% at ";
   printTime();
   for(int i = 0; i < nJobs; i++){
   job(res[i], paramList[i]);
   cout << setprecision(4) << 100*(double)(i+1)/(double)nJobs << " \% at ";
   printTime();
   if(printToFile){
   out->line(res[i]);
   }
   }
   }
   */

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

void MultiThread::setAuxFile(char * argv0, string fileName){
  existsAuxFile = true;
  auxFileName = fileName;
  string path(argv0);
  string dir = path.substr(0, path.find_last_of('/') + 1);
}
