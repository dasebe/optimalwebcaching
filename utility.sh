#!/bin/bash

TIMECMD='/usr/bin/time -f "%e %S %M %C"'

for i in msr_proj_0.tr msr_proj_1.tr msr_proj_2.tr msr_src1_0.tr traceHK100m.tr traceUS100m.tr w100m.tr w500m.tr; do
        echo "UtilityKnapsack/utility ~/LemonMCF/ntrace/${i} > nsolution/sol_utilityknapsack_${i}.log"
done | parallel -j 10 --tmpdir /home/dberger/tmp



