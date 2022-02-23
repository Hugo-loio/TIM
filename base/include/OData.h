#ifndef ODATA_h
#define ODATA_h

#include <string>
#include <fstream>
#include <vector>
using namespace std;

class OData{
  public:
    OData(char * argv0, string fName); //Output argv[0] to argv0
    ~OData();

    void line(string line);
    void data(double ** data, int dim, int nPoints);
    void data(vector<vector<double>> data);

  private:
    ofstream f;
};

#endif
