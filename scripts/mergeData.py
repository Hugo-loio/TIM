import helper as hp 
import numpy as np
import sys
np.set_printoptions(threshold = sys.maxsize)

rawData = []
for fName in sys.argv[1:-1]:
    rawData.append(hp.readfile(fName))

if(len(rawData) == 0):
    print("No input files provided\n")
    quit()

mergedData = rawData[0]
for data in rawData[1:]:
    mergedData = np.concatenate((mergedData, data[1:]), axis = 0)

hp.writeToFile(sys.argv[-1], mergedData)
