#!/bin/bash

set -e

tmpdir=`mktemp -d`

bcftools=/illumina/thirdparty/bcftools/bcftools-1.5/bcftools

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"/.././test/smoke_test_output/

echo running FORMAT/PS behaviour test....

echo  ./test/test2/NA12877_S1.vcf.gz > gvcfs.txt
echo  ./test/test2/NA12878_S1.vcf.gz >> gvcfs.txt
#./bin/gvcfgenotyper -l gvcfs.txt -f ./test/test2/test2.ref.fa | $bcftools query -f '%CHROM %POS[ %PS]\n' -t chr1:70600-70777 > ${tmpdir}/ps_test1.observed
#diff ${tmpdir}/ps_test1.observed ${DIR}/ps_test1.expected
#echo FORMAT/PS behaviour test passed

echo running set_region test
./bin/gvcfgenotyper -l gvcfs.txt -f ./test/test2/test2.ref.fa -r chr1:90000-95000 | $bcftools query -f '%CHROM %POS %REF %ALT\n' > ${tmpdir}/region_test1.observed
diff ${tmpdir}/region_test1.observed ${DIR}/region_test1.expected
echo set_region test passed

echo running set_region test with edge case
ls ./test/test3/NA1287?_S1.vcf.gz > gvcfs.txt
./bin/gvcfgenotyper -l gvcfs.txt -f ./test/test3/test3.ref.fa -r chr1:56680-60000 -f ./test/test3/test3.ref.fa | $bcftools query -f '[%CHROM %POS %REF %ALT %SAMPLE %GT\n]' > ${tmpdir}/region_test2.observed
diff ${tmpdir}/region_test2.observed ${DIR}/region_test2.expected
echo set_region test passed

rm -rf $tmpdir

