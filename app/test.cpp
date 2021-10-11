// basic file operations
#include <iostream>
#include "odata.h"

using namespace std;

int main (int arc, char ** argv) {
  cout << arc << endl;
  cout << argv[0] << endl;

  odata od(argv[0], "test.txt");
  od.line("test the line function");


  double ** dat = new double*[2];
  dat[0] = new double[3];
  dat[0][0] = 1;
  dat[0][1] = 2;
  dat[0][2] = 3;

  dat[1] = new double[3];
  dat[1][0] = 3;
  dat[1][1] = 2;
  dat[1][2] = 4;

  od.data(dat, 2, 3); 

  vector<vector<double>> dat2; 
  vector<double> dat20;
  dat20.push_back(1);
  dat20.push_back(2);
  dat20.push_back(3);
  vector<double> dat21;
  dat21.push_back(3);
  dat21.push_back(2);
  dat21.push_back(1);
  dat2.push_back(dat20);
  dat2.push_back(dat21);

  od.data(dat2);

  return 0;
}
