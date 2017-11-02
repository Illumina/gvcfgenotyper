#!/usr/bin/env bash
#
#

set -o nounset

cppcheck=/illumina/thirdparty/cppcheck/cppcheck-1.69/cppcheck

thisDir=$(dirname $0)
cxx_base_dir=$thisDir/../../

cd $cxx_base_dir 
${cppcheck} -v --enable=all --std=c++11 src/ include/ 2> err.log

