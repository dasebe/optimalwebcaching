 #!/bin/bash

for i in w.tr5m traceHK.tr5m traceUS.tr5m msr_proj_0.tr5m msr_proj_2.tr5m msr_src1_0.tr5m; do #msr_proj_1.tr5m msr_usr_1.tr5m memcachier.tr5m
    for size in $(seq 24 2 42) $(seq 10 2 22); do
        bsize=$((2**${size}))
        echo "./localratio ntrace/${i}_500k ${bsize} | tee nsolution/sol_localratio_${i}_${size}.log"
    done
done | parallel -j 8 --tmpdir tmp

