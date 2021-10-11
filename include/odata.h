#ifndef ODATA_h
#define ODATA_h

#include <string>
#include <fstream>
#include <vector>
using namespace std;

class odata{
  public:
    odata(char * argv0, string fname); //Output argv[0] to argv0
    ~odata();

    void line(string line);
    void data(double ** data, int dim, int npoints);
    void data(vector<vector<double>> data);

  private:
    ofstream f;
};

#endif
