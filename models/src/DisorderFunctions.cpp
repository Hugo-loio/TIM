#include "DisorderFunctions.h"
#include <random>

using namespace std;

double uniform(double from, double to){
  random_device dev;
  mt19937 generator(dev());
  uniform_real_distribution<double> uni(from, to);
  return uni(generator);
}

complex<double> probDisorder2D(Hop & hop, double prob){
  bool isIntraHop = true;
  for(int i = 0; i < 1; i++){
    if(hop.getN()[i] != 0){
      isIntraHop = false;
    }
  }
  if(isIntraHop){
    if(uniform(0,1) < prob){
      return (complex<double>)1/hop.getHop();
    }
    else{
      return hop.getHop();
    }
  }
  else{
    return hop.getHop();
  }
}

complex<double> probDisorder3D(Hop & hop, double prob){
  bool isIntraHop = true;
  for(int i = 0; i < 2; i++){
    if(hop.getN()[i] != 0){
      isIntraHop = false;
    }
  }
  if(isIntraHop){
    if(uniform(0,1) < prob){
      return (complex<double>)1/hop.getHop();
    }
    else{
      return hop.getHop();
    }
  }
  else{
    return hop.getHop();
  }
}

complex<double> sotaiDisorder(Hop & hop, double w){
  bool isIntraHop = true;
  for(int i = 0; i < 2; i++){
    if(hop.getN()[i] != 0){
      isIntraHop = false;
    }
  }
  if(isIntraHop){
    return hop.getHop() + w*uniform(-0.5,0.5);
  }
  else{
    return hop.getHop();
  }
}
