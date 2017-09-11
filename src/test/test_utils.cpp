#include "gtest/gtest.h"
#include "utils.hpp"
#include "common.hpp"


TEST(Utils,util_tests)
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

    ASSERT_TRUE(bcf1_less_than(record1,record1));
    ASSERT_FALSE(bcf1_equal(record1,record2));
    ASSERT_TRUE(bcf1_less_than(record1,record2));
    ASSERT_FALSE(bcf1_greater_than(record1,record2));        
}
