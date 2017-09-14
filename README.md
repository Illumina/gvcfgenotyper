# insert catchy name here
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
./gvcfmerge data/tiny.ref.fa data/NA*.vcf.gz | bcftools view -Ob -o test.bcf
bcftools index test.bcf
```

## run tests:

```
./test_gvcfmerge 
./test_gvcfmerge --gtest_filter=UtilTest.comparators
```
