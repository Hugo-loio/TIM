#include <iostream>
#include "DisorderedSSH.h"
#include "OData.h"
#include "MultiThread.h"
#include "AuxFunctions.h"

void threadPol(int * l, int nSamples, vector<double> & res, vector <double> params){
  DisorderedSSH ssh(params[1]);
  ssh.setSize(l);
  ssh.setW(params[0]);
  res.push_back(params[0]);
  res.push_back(params[1]);

  double pol;
  for(int i = 0; i < nSamples; i++){
    ssh.generateDisorder();
    pol = ssh.polarization();
    res.push_back(pol);
  }
}

void threadIPR(int * l, int nSamples, int nStates, vector<double> & res, vector <double> params){
  DisorderedSSH ssh(params[1]);
  ssh.setSize(l);
  ssh.setW(params[0]);
  res.push_back(params[0]);
  res.push_back(params[1]);

  double ipr;
  for(int i = 0; i < nSamples; i++){
    try{
      ssh.generateDisorder();
      ipr = ssh.ipr(nStates);
      res.push_back(ipr);
    }
    catch(const runtime_error & error){
      cout << "Failed at intra = "  << params[1] << "  w = " << params[0] <<  endl;
      i--;
    }
  }
}

void threadEnt(int * l, int nSamples, vector<double> & res, vector <double> params){
  DisorderedSSH ssh(params[1]);
  ssh.setSize(l);
  ssh.setW(params[0]);
  res.push_back(params[0]);
  res.push_back(params[1]);

  double s;
  int e;
  for(int i = 0; i < nSamples; i++){
    ssh.generateDisorder();
    for(e = 0; e < 3; e++){
      s = ssh.entanglement(e);
      res.push_back(s);
    }
    /*
       s = ssh.entanglement(1);
       res.push_back(s);
       */
  }
}

void pol1(vector<double> & res, vector<double> params){
  int l[1] = {50};
  threadPol(l, 200, res, params);
}

void pol2(vector<double> & res, vector<double> params){
  int l[1] = {100};
  threadPol(l, 200, res, params);
}

void pol3(vector<double> & res, vector<double> params){
  int l[1] = {200};
  threadPol(l, 200, res, params);
}

void pol4(vector<double> & res, vector<double> params){
  int l[1] = {64};
  threadPol(l, 400, res, params);
}

void ipr1(vector<double> & res, vector<double> params){
  int l[1] = {100};
  threadIPR(l, 25, 10, res, params);
}

void ipr2(vector<double> & res, vector<double> params){
  int l[1] = {200};
  threadIPR(l, 25, 10, res, params);
}

void ipr3(vector<double> & res, vector<double> params){
  int l[1] = {500};
  threadIPR(l, 25, 10, res, params);
}

void ent1(vector<double> & res, vector<double> params){
  int l[1] = {50};
  threadEnt(l, 200, res, params);
}

void ent2(vector<double> & res, vector<double> params){
  int l[1] = {100};
  threadEnt(l, 200, res, params);
}

void ent3(vector<double> & res, vector<double> params){
  int l[1] = {200};
  threadEnt(l, 200, res, params);
}

void ent4(vector<double> & res, vector<double> params){
  int l[1] = {20};
  threadEnt(l, 200, res, params);
}

void ent5(vector<double> & res, vector<double> params){
  int l[1] = {64};
  threadEnt(l, 400, res, params);
}

int main (int argc, char ** argv) {
  int threadNumber = 8;
  int version = 0;
  if(argc > 1){
    threadNumber = stoi(argv[1]);
  }
  if(argc > 2){
    version = stoi(argv[2]);
  }

  vector<vector<double>> paramList;
  int nPoints = 30;
  for(int i = 0; i < nPoints; i++){
    for(int e = 0; e < nPoints; e++){
      vector<double> param; 
      param.push_back(6*((double)i+0.5)/(double)nPoints);
      param.push_back(2*((double)e+0.5)/(double)nPoints);
      paramList.push_back(param);
    }
  }

  //run(pol1, paramList, threadNumber, argv[0], "phaseDiagramSSHpol_L50", 0, 0);
  //run(pol2, paramList, threadNumber, argv[0], "phaseDiagramSSHpol_L100", 0, 0);
  //run(pol3, paramList, threadNumber, argv[0], "phaseDiagramSSHpol_L200", 0, 0);
  //run(pol4, paramList, threadNumber, argv[0], "phaseDiagramSSHpol_L64", 0, 0);

  //runSingleThread(ipr1, paramList, argv[0], "phaseDiagramSSHipr_L100_n10", version, 0);
  //runSingleThread(ipr2, paramList, argv[0], "phaseDiagramSSHipr_L200_n10", version, 0);
  //runSingleThread(ipr3, paramList, argv[0], "phaseDiagramSSHipr_L500_n10", version, 0);

  //run(ent1, paramList, threadNumber, argv[0], "phaseDiagramSSHent_L50", 0, 0);
  //run(ent2, paramList, threadNumber, argv[0], "phaseDiagramSSHent_L100", 0, 0);
  //run(ent3, paramList, threadNumber, argv[0], "phaseDiagramSSHent_L200", 0, 0);
  //run(ent4, paramList, threadNumber, argv[0], "phaseDiagramSSHent_L20", 0, 0);
  run(ent5, paramList, threadNumber, argv[0], "phaseDiagramSSHent_L64", 0, 0);
}
