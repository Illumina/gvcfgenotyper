* I am trying to merge a large number of GVCF files and, after opening several files, GVCFgenotyper dies with the error "problem opening ..."

If the file exists and is readable, checkout the max number of file handles that you can open at the same time (ulimit -a).

* How do I create site-only vcf file from the aggregated multi-sample gvcf?

Using bcftools: bcftools view -Ou -G | bcftools norm -m -any -Ou | bcftools view -Oz -o sites.vcf.gz --threads 4
