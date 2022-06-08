#include <iostream>
#include <iomanip>
#include <thread>
#include "DisorderedSOTAI.h"
#include "OData.h"

void scanQuadrupole(int * l, char * argv0, string name){
  int nPoints = 20;
  int nAvg = 10;
  double m = 1.1;

  OData out(argv0, name);
  double delta = 9/(double)nPoints;
  vector<vector<double>> data;
  vector<double> w;
  vector<vector<double>> quad;
  DisorderedSOTAI sotai(m);

  for(int e = 0; e < nAvg; e++){
    quad.push_back(w);
  }

  for(int i = 0; i <= nPoints; i++){
    w.push_back(i*delta);
    sotai.setW(i*delta);
    cout << i << endl;
    for(int e = 0; e < nAvg; e++){
      sotai.generateDisorder();
      quad[e].push_back(sotai.getQuadrupoleManyBody(l));
    }
  }
  data.push_back(w);
  for(int e = 0; e < nAvg; e++){
    data.push_back(quad[e]);
  }
  out.data(data);
}

//Multithreading solution
void threadCalcQuad(int nAvg, int nThreads, int threadIndex, int * l, OData & out){
  int nPoints = 21;
  double m = 1.1;
  int pointsPerThread = nPoints/nThreads;
  int addPoints = 0;
  if(nPoints % nThreads > threadIndex){
    pointsPerThread++;
  }
  else{
    addPoints = nPoints % nThreads;
  }
  int nStart = pointsPerThread*threadIndex + addPoints;
  int nEnd = pointsPerThread*(threadIndex + 1) + addPoints;
  double delta = 9/(double)(nPoints - 1);

  vector<vector<double>> data;
  vector<double> w;
  vector<vector<double>> quad;
  DisorderedSOTAI sotai(m);

  for(int e = 0; e < nAvg; e++){
    quad.push_back(w);
  }

  for(int i = nStart; i < nEnd; i++){
    w.push_back(i*delta);
    sotai.setW(i*delta);
    cout << i << " " << threadIndex << endl;
    for(int e = 0; e < nAvg; e++){
      sotai.generateDisorder();
      quad[e].push_back(sotai.getQuadrupoleManyBody(l));
    }
  }

  data.push_back(w);
  for(int e = 0; e < nAvg; e++){
    data.push_back(quad[e]);
  }
  out.data(data);
}



void scanQuadrupoleMulti(int nAvg, int * l, char * argv0, string name){
  int nThreads = 8;
  vector<thread> thrds;
  OData out(argv0, name);

  for(int i = 0; i < nThreads; i++){
    thrds.push_back(thread(threadCalcQuad, nAvg, nThreads, i, l, ref(out)));
  }

  for(int i = 0; i < nThreads; i++){
    thrds[i].join();
  }

}


int main (int arc, char ** argv) {
  //SOTAI
  int l2D1[2] = {10,10};
  int l2D2[2] = {20,20};
  int l2D3[2] = {30,30};
  cout << "Disordered SOTAI" << endl;
  //scanQuadrupoleMulti(40, l2D1, argv[0], "phaseDiagramSOTAI_10x10_m1.1.dat");
  //scanQuadrupoleMulti(40, l2D2, argv[0], "phaseDiagramSOTAI_20x20_m1.1.dat");
  //scanQuadrupoleMulti(1, l2D3, argv[0], "phaseDiagramSOTAI_30x30_m1.1.dat");
}
