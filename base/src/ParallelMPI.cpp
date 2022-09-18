#include "ParallelMPI.h"
#include <iomanip>

ParallelMPI::ParallelMPI(int * argc, char *** argv){
  MPI_Init(argc, argv);
  MPI_Comm_size(MPI_COMM_WORLD, & nProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, & rank);
}

ParallelMPI::~ParallelMPI(){
  if(areParamsSet){
    deleteParamList();
  }
  if(isFileSet){
    delete out;
  }
  MPI_Finalize();
}

void ParallelMPI::setSamples(int nSamples){
  this->nSamples = nSamples;
}

void ParallelMPI::setFile(char * argv0, string fileName){
  this->fileName = fileName;
  if(isFileSet){
    delete out;
  }
  out = new OData(argv0, fileName);
  isFileSet = true;
}

void ParallelMPI::setParamList(double ** pList, int pSize, int nP){
  if(areParamsSet){
    deleteParamList();
  }
  nParams = nP;
  paramSize = pSize;
  paramList = new double * [nParams];
  int e;
  for(int i = 0; i < nParams; i++){
    paramList[i] = new double[paramSize];
    for(e = 0; e < paramSize; e++){
      paramList[i][e] = pList[i][e];
    }
  }
  areParamsSet = true;
}

void ParallelMPI::setParamList(vector<vector<double>> pList){
  if(areParamsSet){
    deleteParamList();
  }
  nParams = pList.size();
  paramSize = pList[0].size();
  paramList = new double * [nParams];
  int e;
  for(int i = 0; i < nParams; i++){
    paramList[i] = new double[paramSize];
    for(e = 0; e < paramSize; e++){
      paramList[i][e] = pList[i][e];
    }
  }
  areParamsSet = true;
}

void ParallelMPI::deleteParamList(){
  for(int i = 0; i < nParams; i++){
    delete[] paramList[i];
  }
  delete[] paramList;
}

void ParallelMPI::setJob(void (*job) (double * res, double * params), int resSize){
  this->resSize = resSize;
  this->job = job;
  isJobSet = true;
}

void ParallelMPI::run(){
  if(!isJobSet){
    cout << "ERROR: " <<  __PRETTY_FUNCTION__ << " : job not set" << endl;
    return;
  }
  if(!isFileSet){
    cout << "ERROR: " <<  __PRETTY_FUNCTION__ << " : file not set" << endl;
    return;
  }
  if(!areParamsSet){
    cout << "ERROR: " <<  __PRETTY_FUNCTION__ << " : parameters not set" << endl;
    return;
  }

  jobCount = 0;

  tStart = chrono::high_resolution_clock::now();

  if(rank == 0){
    //Master
    fullRes.clear();
    int nJobs = nParams*nSamples;
    cout << "Printing to file: " << fileName << endl;
    cout << "Running " << nJobs << " jobs in " << nProcs - 1 << " processes" << endl;

    int * done = new int[nParams];
    int * sent = new int[nParams];
    int e;
    vector<double> params;
    for(int i = 0; i < nParams; i++){
      done[i] = 0;
      sent[i] = 0;
      params.clear();
      for(e = 0; e < paramSize; e++){
	params.push_back(paramList[i][e]);
      }
      fullRes.push_back(params);
    }
    int countSent = 0;
    int countDone = 0;

    if(nJobs < nProcs - 1){
      nProcs = nJobs + 1;
    }

    MPI_Request req;
    for(int i = 0; i < nProcs - 1; i++){
      e = i/nSamples;
      MPI_Isend(&e, 1, MPI_INT, i + 1, 0, MPI_COMM_WORLD, &req);
      MPI_Request_free(&req);
      sent[e]++;
      if(sent[e] == nSamples){
	countSent++;
      }
    }

    int pIndex, proc;
    double * res = new double[resSize];
    MPI_Status status;
    int flag;
    while(countDone < nParams){
      flag = 0;
      while(1){
	MPI_Iprobe(MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);
	if(flag){
	  break;
	}
	else{
	  sleep(0.3);
	}
      }
      MPI_Recv(&pIndex, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
      proc = status.MPI_SOURCE;
      MPI_Recv(res, resSize, MPI_DOUBLE, proc, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      for(int i = 0; i < resSize; i++){
	fullRes[pIndex].push_back(res[i]);
      }
      done[pIndex]++;
      if(done[pIndex] == nSamples){
	out->clear();
	out->data(fullRes, 0);
	countDone++;
	cout << setprecision(4) << 100*(double)countDone/(double)nParams << " \% at ";
	printTime();
      }
      else if(printEachSamp){
	out->clear();
	out->data(fullRes, 0);
      }
      if(countSent < nParams){
	MPI_Isend(&countSent, 1, MPI_INT, proc, 0, MPI_COMM_WORLD, &req);
	MPI_Request_free(&req);
	sent[countSent]++;
	if(sent[countSent] == nSamples){
	  countSent++;
	}
      }
      else{
	MPI_Isend(&countSent, 1, MPI_INT, proc, 2, MPI_COMM_WORLD, &req);
	MPI_Request_free(&req);
      }
    }
    delete[] sent;
    delete[] done;
    delete[] res;
  }
  else{
    //Slaves
    double * res = new double[resSize];
    int pIndex;
    int countSent;
    int flag;
    int end = 0;
    MPI_Request req;
    while(1){
      flag = 0;
      while(!flag){
	MPI_Iprobe(0, 0, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);
	MPI_Iprobe(0, 2, MPI_COMM_WORLD, &end, MPI_STATUS_IGNORE);
	if(end){
	  delete[] res;
	  cout << "Process " << rank << " did " << jobCount << " jobs" << endl;
	  printTimePerJob();
	  MPI_Recv(&countSent, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  return;
	}
      }
      MPI_Recv(&pIndex, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      job(res, paramList[pIndex]);
      jobCount++;

      MPI_Isend(&pIndex, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &req);
      MPI_Request_free(&req);
      MPI_Isend(res, resSize, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &req);
      MPI_Request_free(&req);
    }
  }
}

void ParallelMPI::printTime(){
  auto tNow = chrono::high_resolution_clock::now();
  auto dH = chrono::duration_cast<chrono::hours>(tNow-tStart);
  auto dM = chrono::duration_cast<chrono::minutes>(tNow-tStart);
  auto dS = chrono::duration_cast<chrono::seconds>(tNow-tStart);
  int h = dH.count();
  int m = dM.count() % 60;
  int s = dS.count() % 60;
  cout << h << ":" << m << ":" << s << endl;
}

void ParallelMPI::printTimePerJob(){
  auto tNow = chrono::high_resolution_clock::now();
  auto dS = chrono::duration_cast<chrono::seconds>(tNow-tStart);
  int s = dS.count() / jobCount;
  int m = (s/60) % 60;
  int h = s/3600;
  s %= 60;
  cout << "Each job in process " << rank << " took an average time of " << h << ":" << m << ":" << s << endl;
}
