#include "gtest/gtest.h"
#include "utils.hpp"
#include "common.hpp"


TEST(UtilTest,comparators)
{
    std::string gvcf_file_name = g_testenv->getBasePath() + "/data/NA12877.tiny.vcf.gz";    
    
    bcf_hdr_t *hdr = bcf_hdr_read(hts_open(gvcf_file_name.c_str(),"r"));
    
    bcf1_t *record1 = bcf_init1();
    record1->rid = 2;
    record1->pos = 92;
    bcf_update_alleles_str(hdr, record1, "TATTAAGATTG,AAGGTTT");
    bcf1_t *record2 = bcf_init1();
    record2->rid = 2;
    record2->pos = 2022;
    bcf_update_alleles_str(hdr, record2, "G,C");

    ASSERT_FALSE(bcf1_less_than(record1,record1));
    ASSERT_TRUE(bcf1_less_than(record1,record2));
    ASSERT_FALSE(bcf1_equal(record1,record2));
    ASSERT_TRUE(bcf1_less_than(record1,record2));
    ASSERT_FALSE(bcf1_greater_than(record1,record2));        

    record2->pos=record1->pos=2399;
    bcf_update_alleles_str(hdr, record1, "C,CTTTTTT");	
    bcf_update_alleles_str(hdr, record2, "CTTTTT,C");
    ASSERT_TRUE(bcf1_less_than(record1,record2));
    ASSERT_FALSE(bcf1_less_than(record2,record1));
}

TEST(UtilTest,GenotypeIndex)
{
    ASSERT_EQ(0,get_gl_index(0,0));
    ASSERT_EQ(1,get_gl_index(0,1));
    ASSERT_EQ(2,get_gl_index(1,1));
    ASSERT_EQ(3,get_gl_index(0,2));
    ASSERT_EQ(4,get_gl_index(1,2));
    ASSERT_EQ(5,get_gl_index(2,2));
}

TEST(UtilTest,phred)
{
    ASSERT_EQ(0,phred(1.0));
    ASSERT_EQ(30,phred(.001));
    ASSERT_FLOAT_EQ(1.0,unphred(0));
    ASSERT_FLOAT_EQ(.001,unphred(30));
}

