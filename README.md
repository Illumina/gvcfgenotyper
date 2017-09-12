# insert catchy name here
merging and genotyping tool for illumina gvcfs

## building

### debug:

```
mkdir debug
cd debug/
cmake ../  -DCMAKE_BUILD_TYPE=Debug
make
```

### release:

```
mkdir debug
cd debug/
cmake ../  -DCMAKE_BUILD_TYPE=Debug
make
```

## Running:

```
./gvcfmerge data/tiny.ref.fa data/NA*.vcf.gz | bcftools view -Ob -o test.bcf
bcftools index test.bcf
```
