
* I am trying to merge a large number of GVCF files and, after opening several files, GVCFgenotyper dies with the error "problem opening ..."

If the file exists and is readable, checkout the max number of file handles that you can open at the same time (ulimit -a).

* How do I create site-only vcf file from the aggregated multi-sample gvcf?

Using bcftools: bcftools view -Ou -G | bcftools norm -m -any -Ou | bcftools view -Oz -o sites.vcf.gz --threads 4

* Eror message: "VCF record did not match the reference at sample..."

This is caused by a bug in an (outdated) version of strelka2. We recommend re-analyzing your samples with a new version of strelka, since you will also benefit from recent improvements in calling accuracy.
If you still want to continue, you can use the undocumented command line option "--ignore-non-matching-ref" to ignore this error.

