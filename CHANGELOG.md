# Changelog
All notable changes to this project will be documented in this file.

# 2018-03-02
- added Fisher's exact test for strand bias (INFO/FS)
- added INFO/GQ_MEDIAN and INFO/GQX_MEDIAN to output
- updated to htslib-1.7
- improved some error messages when non-GVCFs are provided as input

# 2018-02-07
- first public release

## [Unreleased] - 2017-12-20
- replaced CMake with good old Make

## [Unreleased] - 2017-11-01
- updated to htslib-1.6
- first version that produces multi-allelic output

## [Unreleased] - 2017-10-27
- heavily modified mnp_split function such that it will correctly handle (all?) multi-allelic MNPs
- introduced the Genotype class which stores and manipulates our canonical FORMAT fields (currently GQ,DP,DPF,AD,PL)
- buffers now store a fixed number of base pairs, not  number of variants
- added -r argument to set regions
- substantial performance improvements via better buffer memory management
- added handling of hemizygous genotypes
- fixed a bug that caused the buffers to become out of sync across chromosomes
- fixed a small memory leak in GVCFReader closing its bcf_sr_reader
- removed unneccessary vt code

## [0.0.0] - 2017-09-12
- first working release
