# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased] - 2016-11-01
-updated to htslib-1.6
-first version that produces multi-allelic output

## [Unreleased] - 2016-10-27
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
