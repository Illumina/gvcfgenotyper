# Adding additional flags

We recommend bcftools for any post-aggregation of the gvcf. To give an example, bcftools comes with the handy plugin fill-tags to compute various INFO tags such as Minor-Allele-Frequency (MAF) or the perform a test of Hardy-Weinberg equilibrium. 

Usage is straightforward:

```
bcftools +fill-tags in.bcf -Ob -o out.bcf -- -S sample-group.txt -t HWE
```

This performs a test for Hardy-Weinberg equilibrium (HWE).

More info on bcftools [plugins](https://samtools.github.io/bcftools/howtos/plugins.html).

Note that gvcfgenotyper already computes some of the bcftools tags such as AC, AN, DP_MEDIAN and DP_HIST_ALT.

