#include <iostream>
#include <iomanip>
#include "DisorderedSOTAI.h"
#include "OData.h"

void scanQuadrupole(int nPoints, int nAvg, double m, int * l, char * argv0, string name){
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


int main (int arc, char ** argv) {
  //SOTAI
  int l2D1[2] = {10,10};
  int l2D2[2] = {20,20};
  cout << "Disordered SOTAI" << endl;
  //scanQuadrupole(20, 10, 1.1, l2D1, argv[0], "phaseDiagramSOTAI_10x10_m1.1.dat");
  scanQuadrupole(20, 10, 1.1, l2D2, argv[0], "phaseDiagramSOTAI_20x20_m1.1.dat");
}
