#!/bin/bash

TIMECMD='/usr/bin/time -f "%e %S %M %C"'

for i in "traceUS100m.tr"; do
    for size in $(seq 28 1 42); do
        bsize=$((2**${size}))
        for j in 500000 10000 1000000; do
            echo "${TIMECMD} -o ntime/lnslbound_${i}_${j}_${size}.log lnslbound/lbound ../LemonMCF/ntrace/$i ${bsize} 4 $j nsolution/dvars_lnslbound_${i}_${j}_${size}.log | tee nsolution/sol_lnslbound_${i}_${j}_${size}.log"
        done

    done
done | parallel -j 36 --tmpdir /home/dberger/tmp

# =>/usr/bin/time -f "%e %S %M %C" sleep 2
# reals, cput, mrss, cmd
