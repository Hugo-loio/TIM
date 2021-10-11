#include "odata.h"

odata::odata(char * argv0, string fname){
  string path(argv0);
  string dir = path.substr(0,path.find_last_of('/') + 1);
  f.open(dir + "data/" + fname);
}

odata::~odata(){
  f.close();
}

void odata::line(string line){
  f << line << endl;
}

void odata::data(double ** data, int dim, int npoints){
  for(int i = 0; i < npoints; i++){
    for(int e = 0; e < dim; e++){
      f << data[e][i] << " " ; 
    }
    f << endl;
  }
}

void odata::data(vector<vector<double>> data){
  if(data.size() != 0){
    for(int i = 0; i < data[0].size(); i++){
      for(int e = 0; e < data.size(); e++){
	f << data[e][i] << " " ; 
      }
      f << endl;
    }
  }
}
