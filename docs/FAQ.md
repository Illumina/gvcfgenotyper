
* I am trying to merge a large number of GVCF files and, after opening several files, GVCFgenotyper dies with the error "problem opening ..."

If the file exists and is readable, checkout the max number of file handles that you can open at the same time (ulimit -a).

* How do I create site-only vcf file from the aggregated multi-sample gvcf?

Using bcftools: bcftools view -Ou -G | bcftools norm -m -any -Ou | bcftools view -Oz -o sites.vcf.gz --threads 4

* Eror message: "VCF record did not match the reference at sample..."

This is caused by a bug in an (outdated) version of strelka2. We recommend re-analyzing your samples with a newer version of strelka, since you will also benefit from recent improvements in calling accuracy.
If you still want to continue, you can use the undocumented command line option "--ignore-non-matching-ref" to ignore this error.

* How does the varant decomposition in gvcfgenotyper work?

Variants are decomposed and normalized on the fly using code that borrows from "bcftools norm". This process comprises left-alignment and normalization of indels, check if REF alleles match the reference and splitting 
multiallelic sites into multiple rows. Overlapping records are collapsed and sanitized.

* How can I see if a variant PASSed filters in an individual sample?

The current behavior is that the PASS column from the single sample gvcf is propagated into FORMAT/FT in the multisample gvcf. 
Only "." is translated into PASS, all other tags are copied as is.

* How can I get gvcfgenotyper output into Hail?

Check out the instruction in the subdiretory "hail".

