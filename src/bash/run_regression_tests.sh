#!/bin/bash

set -e

tmpdir=`mktemp -d`

echo running regression tests
echo tmpdir $tmpdir

for i in test/regression/*.vcf.gz;
do
    echo Testing $i
    echo $i > test.txt
    bin/gvcfgenotyper -f ${i%vcf.gz}fa -l test.txt | grep -A1000 CHROM > ${tmpdir}/$(basename $i).observed 
    diff ${tmpdir}/$(basename $i).observed ${i}.expected
done

rm test.txt
rm -rf $tmpdir

echo "Regression tests passed"

##this is how you build tests. Use with caution!
#     echo $i > test.txt
#     bin/gvcfgenotyper -f ${i%vcf.gz}fa -l test.txt | grep -A1000 CHROM > ${i}.expected

