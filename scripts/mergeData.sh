#!/bin/sh

if [ $# -ne 1 ]; then
  echo 'Give file name argument'
  exit
fi

scriptDir=$(dirname $0)
dataDir=$scriptDir/../build/data/
parts=$(ls $dataDir | grep -e "$1_p.(1)" | sed -e 's/(1)\.dat//g')

if [ -z "$parts" ]; then
  rawData=$(ls $dataDir | grep $1'(')
  (set -x ; python $scriptDir/mergeData.py $rawData $1'.dat')
else
  for part in $parts
  do
    rawData=$(ls $dataDir | grep $part'(')
    (set -x ; python $scriptDir/mergeData.py $rawData $part'.dat')
  done

  file="$1.dat"
  if [ ! -z "$(ls $dataDir | grep $file)" ]; then
    printf "Overwrite contents of $file ? [y/n]"
    read answer 
    if [ "$answer" != "y" ] && [ "$answer" != "Y" ] ; then
      echo "$file not overwritten"
      exit
    fi
  fi

  > $dataDir$file
  for part in $parts
  do
    cat $dataDir$part'.dat' >> $dataDir$file
  done
fi
