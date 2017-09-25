#include "gtest/gtest.h"
#include "GVCFMerger.hpp"
#include "common.hpp"
#include "StringUtil.hpp"

#include <dirent.h>

TEST(GVCFMerger,platinumGenomeTinyTest)
{
    std::vector<std::string> files;
    std::string test_base = g_testenv->getBasePath() + "/data/test2/";
    DIR* dir;
    struct dirent* ent;
    if ((dir = opendir(test_base.c_str())) != NULL)
    {
        while ((ent = readdir(dir)) != NULL)
        {
            const std::string fname = ent->d_name;
            if (stringutil::endsWith(fname, ".vcf.gz"))
            {
                files.push_back(test_base+fname);
            }
        }
        closedir(dir);
    }
    else
    {
        FAIL() << "Directory of test cases was not found at " << test_base;
    }

    int buffer_size=200;

    std::string ref_file_name = g_testenv->getBasePath() + "/data/test2/test2.ref.fa";    
    std::string output_file_name = std::tmpnam(NULL);
//    std::cerr << "Outputting to " << output_file_name << std::endl;
    GVCFMerger g(files,output_file_name,"z",ref_file_name,buffer_size);
    g.write_vcf();    
}
