#!/bin/bash

export MKL_CBWR=AUTO,STRICT
export MKL_DEBUG_CPU_TYPE=5
export MKL_NUM_THREADS=1
export DR_HOOK=1
export DR_HOOK_OPT=prof
export DR_HOOK_SILENT=1


for opt in "" "-use_trans1" 
do

~/SAVE/mpiauto/mpiauto \
 --wrap --wrap-stdeo -np 4 -openmp 1 -- \
 ./bin/ectrans-benchmark-dp \
    -n 16 --nproma 32 \
    -t 79 -f 8 -l 15 -g N64 $opt

if [ "x$opt" != "x" ]
then
  \rm -rf T
  mkdir T
  \mv stdeo.* trans.*.dat T/
else
  \rm -rf F
  mkdir F
  \mv stdeo.* trans.*.dat F/
fi

done

echo
echo

echo "==> DIFF <=="

for f in T/trans.*
do
  b=$(basename $f)
  echo "==> $b <=="
  diff T/$b F/$b
done
