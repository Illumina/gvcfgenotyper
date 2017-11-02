# gvcfgenotyper
merging and genotyping tool for illumina gvcfs


### build:

```
git clone git@git.illumina.com:Bioinformatics/gvcfgenotyper.git
mkdir release
cd release/
cmake ../
make
```

### Running:

This takes about 30 minutes:

```
find /illumina/build/platinumgenomes/builds/hg19/pg_ns6/ -name '*genome.vcf.gz' > gvcfs.txt
ref=/illumina/development/Isis/Genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
time ./gvcfgenotyper -f $ref -l gvcfs.txt -Ob -o pg.ns6.bcf
```

or with some trivial parallelism:

```
for i in {1..22} X;do echo -f $ref -l gvcfs.txt -Ob -o pg.ns6.chr${i}.bcf;done | xargs -l -P 23 ./gvcfgenotyper
```


### debug:

```
mkdir debug
cd debug/
cmake ../  -DCMAKE_BUILD_TYPE=Debug  -DCMAKE_CXX_COMPILER=`which g++` -DCMAKE_C_COMPILER=`which gcc`
make
```

### run tests:

```
##run all tests
make test

##run google tests
./test_gvcfgenotyper

##run smoke tests
bash -e bash/smoke_tests.sh
```
