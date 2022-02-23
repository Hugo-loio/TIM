#include "OData.h"

OData::OData(char * argv0, string fName){
  string path(argv0);
  string dir = path.substr(0,path.find_last_of('/') + 1);
  f.open(dir + "data/" + fName);
}

OData::~OData(){
  f.close();
}

void OData::line(string line){
  f << line << endl;
}

void OData::data(double ** data, int dim, int nPoints){
  for(int i = 0; i < nPoints; i++){
    for(int e = 0; e < dim; e++){
      f << data[e][i] << " " ; 
    }
    f << endl;
  }
}

void OData::data(vector<vector<double>> data){
  if(data.size() != 0){
    for(int i = 0; i < data[0].size(); i++){
      for(int e = 0; e < data.size(); e++){
	f << data[e][i] << " " ; 
      }
      f << endl;
    }
  }
}
