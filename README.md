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

This takes about 30 minutes:

```
find /illumina/build/platinumgenomes/builds/hg19/pg_ns6/ -name '*genome.vcf.gz' > gvcfs.txt
ref=/illumina/development/Isis/Genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
time ./gvcfgenotyper -f $ref -l gvcfs.txt -Ob -o pg.ns6.bcf
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
