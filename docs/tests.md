```
 bcftools view -i 'N_ALT>0' data/NA12877.tiny.vcf.gz | bcftools norm -m -any  | bcftools norm -f data/tiny.ref.fa | vt decompose_blocksub   - | vt sort - | vt uniq - | bcftools query  -f '%POS %REF %ALT\n' | sort -k 1,1n -k 2,2 -k 3,3 | awk -v OFS=":" '{print "chr3",$1,$2,$3}'  > data/NA12877.tiny.vcf.gz.expected 
```
