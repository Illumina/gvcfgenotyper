#include "gtest/gtest.h"
#include "utils.hpp"
#include "common.hpp"
#include "test_utils.h"

TEST(UtilTest, comparators)
{
    std::string gvcf_file_name = g_testenv->getBasePath() + "/data/NA12877.tiny.vcf.gz";

    bcf_hdr_t *hdr = bcf_hdr_read(hts_open(gvcf_file_name.c_str(), "r"));

    bcf1_t *record1 = generate_record(hdr, 2,92, "TATTAAGATTG,AAGGTTT");
    bcf1_t *record2 = generate_record(hdr,2,2022,"G,C");


    ASSERT_TRUE(is_snp(record2));
    ASSERT_FALSE(bcf1_less_than(record1, record1));
    ASSERT_TRUE(bcf1_less_than(record1, record2));
    ASSERT_FALSE(bcf1_equal(record1, record2));
    ASSERT_TRUE(bcf1_less_than(record1, record2));
    ASSERT_FALSE(bcf1_greater_than(record1, record2));
    ASSERT_TRUE(bcf1_leq(record2,record2));


    update_record(hdr,2,2400,"C,CTTTTTT",record1);
    ASSERT_FALSE(is_deletion(record1));
    ASSERT_TRUE(is_insertion(record1));
    ASSERT_FALSE(is_snp(record1));

    update_record(hdr,2,2400,"CTTTTTT,C",record2);
    ASSERT_TRUE(is_deletion(record2));
    ASSERT_FALSE(is_insertion(record2));
    ASSERT_FALSE(is_snp(record2));

    ASSERT_TRUE(bcf1_less_than(record1, record2));
    ASSERT_FALSE(bcf1_less_than(record2, record1));

    update_record(hdr,2,2400,"C,G",record1);
    update_record(hdr,2,2400,"CTTTTT,C",record2);
    ASSERT_TRUE(bcf1_less_than(record1, record2));
    ASSERT_FALSE(bcf1_less_than(record2, record1));

    print_variant(hdr,record1);
    std::cerr << bcf_get_variant_type(record1,1)<<std::endl;
    ASSERT_TRUE(is_snp(record1));
    ASSERT_FALSE(is_insertion(record1));
    ASSERT_FALSE(is_deletion(record1));

    update_record(hdr,2,2400,"C,CGGGG",record1);
    ASSERT_FALSE(is_snp(record1));
    ASSERT_TRUE(is_insertion(record1));
    ASSERT_FALSE(is_deletion(record1));
    ASSERT_FALSE(is_complex(record1));

    update_record(hdr,2,2400,"CTTTTTT,C",record1);
    ASSERT_FALSE(is_snp(record1));
    ASSERT_FALSE(is_insertion(record1));
    ASSERT_TRUE(is_deletion(record1));
    ASSERT_FALSE(is_complex(record1));

    update_record(hdr,2,2400,"CTTTTTT,G",record1);
    ASSERT_FALSE(is_snp(record1));
    ASSERT_FALSE(is_insertion(record1));
    ASSERT_FALSE(is_deletion(record1));
    ASSERT_TRUE(is_complex(record1));

    update_record(hdr,2,2400,"T,C",record1);
    update_record(hdr,2,2400,"T,C,X",record2);
    ASSERT_TRUE(bcf1_leq(record1,record2));

    update_record(hdr,2,2400,"C,CA",record1);
    ASSERT_TRUE(is_insertion(record1));
    ASSERT_FALSE(is_deletion(record1));

    update_record(hdr,2,2400,"CAAAAAAA,CA,C,CAA",record2);
    ASSERT_TRUE(is_deletion(record2));
    ASSERT_FALSE(is_insertion(record2));
    ASSERT_TRUE(bcf1_leq(record1,record2));

    update_record(hdr,2,2400,"G,GGTGTGT",record1);
    update_record(hdr,2,2400,"GGGGTGTGTGT,GGT,G",record2);
    ASSERT_EQ(get_variant_rank(record1),1);
    ASSERT_EQ(get_variant_rank(record2),2);

    ASSERT_TRUE(bcf1_greater_than(record2,record1));
    ASSERT_TRUE(bcf1_less_than(record1,record2));
    ASSERT_TRUE(is_insertion(record1));
    print_variant(record1);
    print_variant(record2);

    update_record(hdr,0,11017,"C,T",record1);
    update_record(hdr,1,10000,"T,A",record2);
    ASSERT_FALSE(bcf1_equal(record2,record1));
    ASSERT_FALSE(bcf1_less_than(record2,record1));
    ASSERT_TRUE(bcf1_greater_than(record2,record1));
    bcf_destroy(record1);
    bcf_destroy(record2);
}

TEST(UtilTest, GenotypeIndex)
{
    ASSERT_EQ(0, get_gl_index(0, 0));
    ASSERT_EQ(1, get_gl_index(0, 1));
    ASSERT_EQ(2, get_gl_index(1, 1));
    ASSERT_EQ(3, get_gl_index(0, 2));
    ASSERT_EQ(4, get_gl_index(1, 2));
    ASSERT_EQ(5, get_gl_index(2, 2));
}

TEST(UtilTest, phred)
{
    ASSERT_EQ(0, phred(1.0));
    ASSERT_EQ(30, phred(.001));
    ASSERT_FLOAT_EQ(1.0, unphred(0));
    ASSERT_FLOAT_EQ(.001, unphred(30));
}

