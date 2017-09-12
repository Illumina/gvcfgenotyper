#include "gtest/gtest.h"
#include "GVCFMerger.hpp"
#include "common.hpp"

TEST(GVCFMerger,readThreeGVCFs)
{
    int buffer_size=200;
    std::vector<std::string> files;
    files.push_back(g_testenv->getBasePath() + "/data/NA12877.tiny.vcf.gz");
    std::string ref_file_name = g_testenv->getBasePath() + "/data/tiny.ref.fa";    
    files.push_back(g_testenv->getBasePath() + "/data/NA12889.tiny.vcf.gz");
    files.push_back(g_testenv->getBasePath() + "/data/NA12890.tiny.vcf.gz");
    GVCFMerger g(files,"/tmp/test.vcf.gz","z",ref_file_name,buffer_size);
    g.write_vcf();    
}
