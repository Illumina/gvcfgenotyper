# gvcfgenotyper
merging and genotyping tool for illumina gvcfs

## building

### debug:

```
mkdir debug
cd debug/
cmake ../  -DCMAKE_BUILD_TYPE=Debug  -DCMAKE_CXX_COMPILER=`which g++` -DCMAKE_C_COMPILER=`which gcc`
make
```

### release:

```
mkdir release
cd release/
cmake ../ -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=`which g++` -DCMAKE_C_COMPILER=`which gcc`
make
```

## Running:

```
$ find /illumina/scratch/kimura/gvcfs/6.19.1.403/ -name '*vcf.gz' > full_gvcfs.txt
$ time ./gvcfmerge -f $b37 -l full_gvcfs.txt -Ob -o test.bcf
Input GVCFs:
/illumina/scratch/kimura/gvcfs/6.19.1.403/NA12891_S1_S1.genome.vcf.gz
/illumina/scratch/kimura/gvcfs/6.19.1.403/NA12885_S1_S1.genome.vcf.gz
/illumina/scratch/kimura/gvcfs/6.19.1.403/NA12889_S1_S1.genome.vcf.gz
/illumina/scratch/kimura/gvcfs/6.19.1.403/NA12877_S1_S1.genome.vcf.gz
/illumina/scratch/kimura/gvcfs/6.19.1.403/NA12884_S1_S1.genome.vcf.gz
/illumina/scratch/kimura/gvcfs/6.19.1.403/NA12893_S1_S1.genome.vcf.gz
/illumina/scratch/kimura/gvcfs/6.19.1.403/NA12881_S1_S1.genome.vcf.gz
/illumina/scratch/kimura/gvcfs/6.19.1.403/NA12879_S1_S1.genome.vcf.gz
/illumina/scratch/kimura/gvcfs/6.19.1.403/NA12886_S1_S1.genome.vcf.gz
/illumina/scratch/kimura/gvcfs/6.19.1.403/NA12892_S1_S1.genome.vcf.gz
/illumina/scratch/kimura/gvcfs/6.19.1.403/NA12887_S1_S1.genome.vcf.gz
/illumina/scratch/kimura/gvcfs/6.19.1.403/NA12888_S1_S1.genome.vcf.gz
/illumina/scratch/kimura/gvcfs/6.19.1.403/NA12878_S1_S1.genome.vcf.gz
/illumina/scratch/kimura/gvcfs/6.19.1.403/NA12883_S1_S1.genome.vcf.gz
/illumina/scratch/kimura/gvcfs/6.19.1.403/NA12880_S1_S1.genome.vcf.gz
/illumina/scratch/kimura/gvcfs/6.19.1.403/NA12882_S1_S1.genome.vcf.gz
/illumina/scratch/kimura/gvcfs/6.19.1.403/NA12890_S1_S1.genome.vcf.gz
Wrote 10035326 variants

real	33m25.025s
user	30m32.536s
sys	2m32.395s
```

## run tests:

```
##run all tests
make test

##run google tests
./test_gvcfmerge

##run smoke tests
bash -e bash/smoke_tests.sh
```
