# gvcfgenotyper 

A utility for merging and genotyping Illumina-style GVCFs. This software is not commercially supported.

This source code is provided under the [Apache License 2.0](https://choosealicense.com/licenses/apache-2.0/#). Copyright (c) 2018, Illumina, Inc. All rights reserved.


This tool provides basic genome VCF (GVCF) merging and genotyping functionality to provide a multisample BCF/VCF suitable for cohort analysis. Variants are normalised and decomposed on-the-fly before merging. Samples that do not have a particular variant have their homozygous reference confidence estimated from the GVCF depth blocks using some simple heuristics.

#### Caution:
This software is in very early development, it is largely functional but may contain bugs.

There are various flavours of GVCF in the wild, this tool only works with the format [produced by Illumina pipelines](https://sites.google.com/site/gvcftools/home/about-gvcf).


### Installation

```
git clone git@git.illumina.com:Bioinformatics/gvcfgenotyper.git
cd gvcfgenotyper/
make
bin/gvcfgenotyper
```

### Running

```
find directory/ -name '*genome.vcf.gz' > gvcfs.txt
time ./gvcfgenotyper -f genome.fa -l gvcfs.txt -Ob -o output.bcf
```

or with some trivial parallelism:

```
for i in {1..22} X;
do 
    echo -f genome.fa -l gvcfs.txt -Ob -o output.chr${i}.bcf;
done | xargs -l -P 23 ./gvcfgenotyper
```

### Known issues

We are working on multi-threading to improve performance.

Homozygous reference confidence works well for SNPs but is less reliable for indels.

Complex variants can occasionally contain primitive alleles called in other samples. We are investigating decomposition approaches for this problem.

### Acknowledgements

This tool depends on [htslib](www.htslib.org) and [googletest](https://github.com/google/googletest). We also borrowed some variant normalisation code from [BCFtools](https://samtools.github.io/bcftools/bcftools.html).
