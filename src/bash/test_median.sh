#tests that INFO/GQ_MEDIAN and INFO/GQX_MEDIAN are being calculated correctly 
 
 find test/test2/*.vcf.gz > list.txt
 bin/gvcfgenotyper -l list.txt -f test/test2/test2.ref.fa | bcftools norm -m -any | bcftools query -f '%CHROM:%POS:%REF:%ALT %GQX_MEDIAN %GQ_MEDIAN\n' > /tmp/median.txt
 diff /tmp/median.txt test/test2/median.expected.txt