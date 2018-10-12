# Using Hail to analyze GVCFgenotyper output

[Hail](https://github.com/hail-is/hail) is an open-soure framework to analyze genomic data. Hail is based on [Spark](https://spark.apache.org/).

Hail allows scientists to perform statistical analyses on multi-sample VCF files, both exploratory (PCA) and model fitting (linear regression, GWAS).

Hail interfaces with popular python libraries such as [scikit-learn](https://spark.apache.org/).

There are some subtle differences between Illumina-style VCF files and the VCF format that Hail expects. We provide a plugins for bcftools that help with
the conversion.

This assumes that you have downloaded and compiled a recent version of [bcftools](https://samtools.github.io/bcftools/):
 
```
bcftools_dir=/path/to/bcftools
cp ilmn2hail.c $bcftools_dir/plugins
pushd $bcftools_dir
make
popd
$bcftools_dir/bin/bcftools view multi.sample.bcf | bcftools +ilmn2hail | bcftools view -Oz -o multi.sample.for.hail.vcf.bgz
```

The [Hail manual](https://hail.is/docs/devel/methods/impex.html#hail.methods.import_vcf) recommends to import block-compressed VCF. 

If you are looking for an example data set to try this out, have a look at the [Polaris](https://github.com/Illumina/Polaris) cohort.

We provide a [jupyter](./import_vcf_hail_0.2.ipynb) notebook with some Hail commands showing an exploratory analysis of variant calls on chromosome 20 of Polaris.

 
