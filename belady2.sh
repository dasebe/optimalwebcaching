#!/bin/bash

TIMECMD='/usr/bin/time -f "%e %S %M %C"'

for j in 100 500; do 
    for i in "memcachier100m.tr" ; do #"w100m.tr" "traceUS100m.tr" "traceHK100m.tr" ; do  #"w100m.tr" "traceUS100m.tr" "traceHK100m.tr"  "w500m.tr"
        for size in $(seq 20 1 29); do #$(seq 30 2 38) $(seq 27 2 41) 28 40
            bsize=$((2**${size}))
            echo "${TIMECMD} -o ntime/belady2_${i}_${j}_${size}.log Belady2/belady2 ~/LemonMCF/ntrace/$i ${bsize} $j > nsolution/sol_belady2_${i}_${j}_${size}.log"
            echo "${TIMECMD} -o ntime/belady2size_${i}_${j}_${size}.log Belady2Size/belady2size ~/LemonMCF/ntrace/$i ${bsize} $j > nsolution/sol_belady2size_${i}_${j}_${size}.log"
            echo "${TIMECMD} -o ntime/belady2sizefrequency_${i}_${j}_${size}.log Belady2SizeFrequency/belady2sizefrequency ~/LemonMCF/ntrace/$i ${bsize} $j > nsolution/sol_belady2sizefrequency_${i}_${j}_${size}.log"
        done
    done
done | parallel -j 3 --tmpdir /home/dberger/tmp

exit

for j in 100 500 1001 2001 4001 8001 160001; do 
    for i in "memcachier100m.tr" ; do #"w100m.tr" "traceUS100m.tr" "traceHK100m.tr" ; do  #"w100m.tr" "traceUS100m.tr" "traceHK100m.tr"  "w500m.tr"
        for size in $(seq 30 1 30); do #$(seq 30 2 38) $(seq 27 2 41) 28 40
            bsize=$((2**${size}))
#            echo "${TIMECMD} -o ntime/belady2_${i}_${j}_${size}.log Belady2/belady2 ~/LemonMCF/ntrace/$i ${bsize} $j > nsolution/sol_belady2_${i}_${j}_${size}.log"
            echo "${TIMECMD} -o ntime/belady2size_${i}_${j}_${size}.log Belady2Size/belady2size ~/LemonMCF/ntrace/$i ${bsize} $j > nsolution/sol_belady2size_${i}_${j}_${size}.log"
#            echo "${TIMECMD} -o ntime/belady2sizefrequency_${i}_${j}_${size}.log Belady2SizeFrequency/belady2sizefrequency ~/LemonMCF/ntrace/$i ${bsize} $j > nsolution/sol_belady2sizefrequency_${i}_${j}_${size}.log"
        done
    done
done | parallel -j 4 --tmpdir /home/dberger/tmp
