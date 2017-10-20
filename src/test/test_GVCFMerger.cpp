#include "gtest/gtest.h"
#include "GVCFMerger.hpp"
#include "common.hpp"
#include "StringUtil.hpp"

#include <dirent.h>

static bcf1_t *generate_record(bcf_hdr_t *hdr,int rid,int pos,const string & alleles)
{
    bcf1_t *ret = bcf_init1();
    ret->rid = rid;
    ret->pos = pos;
    bcf_update_alleles_str(hdr, ret, alleles.c_str());
    return(ret);
}

static bcf_hdr_t *get_header()
{
    std::string gvcf_file_name = g_testenv->getBasePath() + "/data/NA12877.tiny.vcf.gz";
    return(bcf_hdr_read(hts_open(gvcf_file_name.c_str(), "r")));
}

TEST(multiAllele,test1)
{
    int rid=1;
    int pos=99;
    auto hdr = get_header();
    auto rec1 = generate_record(hdr,rid,pos,"C,G");
    auto rec2 = generate_record(hdr,rid,pos,"C,A");
    auto rec3 = generate_record(hdr,rid,200,"CTG,C");
    auto rec4 = generate_record(hdr,rid,pos,"CTGG,C");
    auto rec5 = generate_record(hdr,rid,pos,"C,CAAAAAAAA");

    multiAllele m(rid,pos,hdr);
    ASSERT_EQ(m.allele(rec1),1);
    ASSERT_EQ(m.allele(rec1),1);
    ASSERT_EQ(m.allele(rec2),2);
    ASSERT_EQ(m.allele(rec1),1);
    ASSERT_EQ(m.allele(rec3),0);
    ASSERT_EQ(m.allele(rec2),2);
    ASSERT_EQ(m.allele(rec4),3);
    ASSERT_EQ(m.allele(rec5),4);

    bcf1_t *v = bcf_init1();
    m.collapse(v);
    print_variant(hdr,v);

    auto truth = generate_record(hdr,rid,pos,"CTGG,GTGG,ATGG,C,CAAAAAAAATGG");     //chr1:100:CTGG:GTGG,ATGG,C,CAAAAAAAATGG
    ASSERT_TRUE(bcf1_equal(truth,v));
    bcf_destroy(rec1);
    bcf_destroy(rec2);
    bcf_destroy(rec3);
    bcf_destroy(rec4);
    bcf_destroy(rec5);
    bcf_destroy(v);
}

TEST(GVCFMerger, platinumGenomeTinyTest)
{
    std::vector<std::string> files;
    std::string test_base = g_testenv->getBasePath() + "/data/test2/";
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir(test_base.c_str())) != NULL)
    {
        while ((ent = readdir(dir)) != NULL)
        {
            const std::string fname = ent->d_name;
            if (stringutil::endsWith(fname, ".vcf.gz"))
            {
                files.push_back(test_base + fname);
            }
        }
        closedir(dir);
    }
    else
    {
        FAIL() << "Directory of test cases was not found at " << test_base;
    }

    int buffer_size = 200;

    std::string ref_file_name = g_testenv->getBasePath() + "/data/test2/test2.ref.fa";
    std::string output_file_name = std::tmpnam(NULL);
//    std::cerr << "Outputting to " << output_file_name << std::endl;
    GVCFMerger g(files, output_file_name, "z", ref_file_name, buffer_size);
    g.write_vcf();
}
