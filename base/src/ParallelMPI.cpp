#include "ParallelMPI.h"
#include <iomanip>

ParallelMPI::ParallelMPI(int * argc, char *** argv){
  MPI_Init(argc, argv);
  MPI_Comm_size(MPI_COMM_WORLD, & nProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, & rank);
}

ParallelMPI::~ParallelMPI(){
  if(areParamsSet){
    for(int i = 0; i < nParams; i++){
      delete[] paramList[i];
    }
    delete[] paramList;
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
  out = new OData(argv0, fileName);
  isFileSet = true;
}

void ParallelMPI::setParamList(double ** pList, int pSize, int nP){
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

  tStart = chrono::high_resolution_clock::now();

  if(rank == 0){
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

    int nJobs = nParams*nSamples;
    if(nJobs < nProcs){
      nProcs = nJobs;
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

    int flag, pIndex, proc;
    double * res = new double[resSize];
    MPI_Status status;
    while(countDone < nParams){
      MPI_Recv(&pIndex, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
      proc = status.MPI_SOURCE;
      MPI_Recv(&res, resSize, MPI_DOUBLE, proc, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      for(int i = 0; i < resSize; i++){
	fullRes[pIndex].push_back(res[i]);
      }
      done[pIndex]++;
      if(done[pIndex] == nSamples){
	out->line(fullRes[pIndex]);
	countDone++;
	cout << setprecision(4) << 100*(double)countDone/(double)nParams << " \% at ";
	printTime();
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
    cout << "done: " << countDone << " sent " << countSent << " nParams " << nParams << endl;
    delete[] sent;
    delete[] done;
    delete[] res;
  }
  else{
    double * res = new double[resSize];
    int pIndex;
    int flag;
    int end = 0;
    while(1){
      flag = 0;
      while(!flag){
	MPI_Iprobe(0, 0, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);
	MPI_Iprobe(0, 2, MPI_COMM_WORLD, &end, MPI_STATUS_IGNORE);
	if(end){
	  delete[] res;
	  cout << "Process " << rank << " did " << jobCount << " jobs" << endl;
	  printTimePerJob();
	  return;
	}
      }
      MPI_Recv(&pIndex, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      job(res, paramList[pIndex]);
      jobCount++;

      MPI_Send(&pIndex, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
      MPI_Send(res, resSize, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    }
  }

}

void ParallelMPI::test(){
  if(rank == 0){
    int size = 4*nProcs;
    int sent = 0;
    int done = 0; 
    MPI_Request req;
    int * array = new int[size];
    for(int i = 0; i < size; i++){
      array[i] = i;
    }

    for(int i = 0; i < nProcs - 1; i++){
      MPI_Isend(&array[i], 1, MPI_INT, i + 1, 0, MPI_COMM_WORLD, &req);
      MPI_Request_free(&req);
      sent++;
    }

    int flag, num, proc;
    MPI_Status status;
    while(done < size){
      MPI_Recv(&num, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
      proc = status.MPI_SOURCE;
      cout << "From process: " << proc << " got num " << num << endl;
      done++;
      if(sent < size){
	MPI_Isend(&array[sent], 1, MPI_INT, proc, 0, MPI_COMM_WORLD, &req);
	MPI_Request_free(&req);
	sent++;
      }
      else{
	MPI_Isend(&sent, 1, MPI_INT, proc, 2, MPI_COMM_WORLD, &req);
	MPI_Request_free(&req);
      }
    }
    cout << "done: " << done << " sent " << sent << " size " << size << endl;
  }
  else{
    int num; 
    int flag;
    int end = 0;
    while(1){
      flag = 0;
      while(!flag){
	MPI_Iprobe(0, 0, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);
	MPI_Iprobe(0, 2, MPI_COMM_WORLD, &end, MPI_STATUS_IGNORE);
	if(end){
	  cout << "Process " << rank << " did " << jobCount << " calculations" << endl;
	  return;
	}
      }
      MPI_Recv(&num, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      num*=2;
      jobCount++;
      MPI_Send(&num, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
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
