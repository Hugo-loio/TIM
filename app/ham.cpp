#include <iostream>
#include "OData.h"
#include "DisorderedSOTAI.h"

using namespace std;
using namespace arma;

int main(int argc, char ** argv){
  DisorderedSOTAI sotai(1.1);
  int l[2] = {4,4};
  sotai.setSize(l);
  sotai.setW(0);
  sotai.generateDisorder();
  sotai.printHam(argv[0], "hamSotai_10x10_w0");
}
