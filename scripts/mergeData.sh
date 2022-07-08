#!/bin/sh

if [ $# -ne 1 ]; then
  echo 'Give file name argument'
  exit
fi

scriptDir=$(dirname $0)
dataDir=$scriptDir/../build/data/
rawData=$(ls $dataDir | grep $1'(')

(set -x ; python $scriptDir/mergeData.py $rawData $1'.dat')

