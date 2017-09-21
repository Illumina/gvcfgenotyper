#!/bin/bash

tmpdir=`mktemp -d`

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"/../data/smoke_test_output/

echo running FORMAT/PS behaviour test....

echo  ../data/test2/NA12877_S1.vcf.gz > gvcfs.txt
echo  ../data/test2/NA12878_S1.vcf.gz >> gvcfs.txt
./gvcfmerge -l gvcfs.txt -f ../data/test2/test2.ref.fa | bcftools query -f '%CHROM %POS[ %PS]\n' -t chr1:70600-70777 > ${tmpdir}/ps_test1.observed
diff ${tmpdir}/ps_test1.observed ${DIR}/ps_test1.expected
echo FORMAT/PS behaviour test passed

rm -rf $tmpdir

