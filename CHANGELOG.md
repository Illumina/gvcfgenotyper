# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]
- added handling of hemizygous genotypes
- fixed a bug that caused the buffers to become out of sync across chromosomes
- fixed a small memory leak in GVCFReader closing its bcf_sr_reader
- removed unneccessary vt code

## [0.0.0] - 2017-09-12
- first working release
