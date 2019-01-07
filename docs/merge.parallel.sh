#/bin/bash

ref=/path/to/your/ref/genome.fa
gvcfs=gvcf.list
bin=path/to/gvcfgenotyper

for i in {1..22} X;
do
    echo -r $i -f $ref -l $gvcfs -Ob -o output.chr${i}.bcf;
 done | xargs -l -P 23 $bin

