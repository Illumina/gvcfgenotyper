#tests that FORMAT/FT is identical to single sample FILTER

echo test/test2/NA12877_S1.vcf.gz > list.txt
bin/gvcfgenotyper -l na12877.txt -f test/test2/test2.ref.fa | bcftools query -f '%CHROM\t%POS\t%REF[\t%FT]\n' > 1.txt
bcftools query -f '%CHROM\t%POS\t%REF[\t%FT]\n'  -i  'N_ALT>0' test/test2/NA12877_S1.vcf.gz  > 2.txt
diff 1.txt 2.txt
echo "FORMAT/FT test passed!"
