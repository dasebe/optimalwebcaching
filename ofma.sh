 #!/bin/bash

for i in w.tr5m traceHK.tr5m traceUS.tr5m msr_proj_0.tr5m msr_proj_1.tr5m msr_proj_2.tr5m msr_src1_0.tr5m msr_usr_1.tr5m memcachier.tr5m; do
    for size in $(seq 22 2 42); do
        bsize=$((2**${size}))
        echo "./ofma ntrace/$i ${bsize} | tee nsolution/sol_ofma_${i}_${size}.log"
    done
done | parallel -j 24 --tmpdir tmp

