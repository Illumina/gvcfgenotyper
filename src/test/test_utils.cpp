#include <htslib/vcf.h>
#include "test_helpers.hh"

#include "ggutils.hh"

TEST(UtilTest, comparators)
{
    std::string gvcf_file_name = g_testenv->getBasePath() + "/data/NA12877.tiny.vcf.gz";

    bcf_hdr_t *hdr = bcf_hdr_read(hts_open(gvcf_file_name.c_str(), "r"));

    bcf1_t *record1 = generate_record(hdr, 2,92, "TATTAAGATTG,AAGGTTT");
    bcf1_t *record2 = generate_record(hdr,2,2022,"G,C");


    ASSERT_TRUE(ggutils::is_snp(record2));
    ASSERT_FALSE( ggutils::bcf1_less_than(record1, record1));
    ASSERT_TRUE( ggutils::bcf1_less_than(record1, record2));
    ASSERT_FALSE( ggutils::bcf1_equal(record1, record2));
    ASSERT_TRUE( ggutils::bcf1_less_than(record1, record2));
    ASSERT_FALSE( ggutils::bcf1_greater_than(record1, record2));
    ASSERT_TRUE( ggutils::bcf1_leq(record2,record2));


    update_record(hdr,2,2400,"C,CTTTTTT",record1);
    ASSERT_FALSE(ggutils::is_deletion(record1));
    ASSERT_TRUE(ggutils::is_insertion(record1));
    ASSERT_FALSE(ggutils::is_snp(record1));

    update_record(hdr,2,2400,"CTTTTTT,C",record2);
    ASSERT_TRUE(ggutils::is_deletion(record2));
    ASSERT_FALSE(ggutils::is_insertion(record2));
    ASSERT_FALSE(ggutils::is_snp(record2));

    ASSERT_TRUE( ggutils::bcf1_less_than(record1, record2));
    ASSERT_FALSE( ggutils::bcf1_less_than(record2, record1));

    update_record(hdr,2,2400,"C,G",record1);
    update_record(hdr,2,2400,"CTTTTT,C",record2);
    ASSERT_TRUE( ggutils::bcf1_less_than(record1, record2));
    ASSERT_FALSE( ggutils::bcf1_less_than(record2, record1));

    ggutils::print_variant(hdr,record1);
    std::cerr << bcf_get_variant_type(record1,1)<<std::endl;
    ASSERT_TRUE(ggutils::is_snp(record1));
    ASSERT_FALSE(ggutils::is_insertion(record1));
    ASSERT_FALSE(ggutils::is_deletion(record1));

    update_record(hdr,2,2400,"C,CGGGG",record1);
    ASSERT_FALSE(ggutils::is_snp(record1));
    ASSERT_TRUE(ggutils::is_insertion(record1));
    ASSERT_FALSE(ggutils::is_deletion(record1));
    ASSERT_FALSE(ggutils::is_complex(record1));

    update_record(hdr,2,2400,"CTTTTTT,C",record1);
    ASSERT_FALSE(ggutils::is_snp(record1));
    ASSERT_FALSE(ggutils::is_insertion(record1));
    ASSERT_TRUE(ggutils::is_deletion(record1));
    ASSERT_FALSE(ggutils::is_complex(record1));

    update_record(hdr,2,2400,"CTTTTTT,G",record1);
    ASSERT_FALSE(ggutils::is_snp(record1));
    ASSERT_FALSE(ggutils::is_insertion(record1));
    ASSERT_FALSE(ggutils::is_deletion(record1));
    ASSERT_TRUE(ggutils::is_complex(record1));

    update_record(hdr,2,2400,"T,C",record1);
    update_record(hdr,2,2400,"T,C,X",record2);
    ASSERT_TRUE( ggutils::bcf1_leq(record1,record2));

    update_record(hdr,2,2400,"C,CA",record1);
    ASSERT_TRUE(ggutils::is_insertion(record1));
    ASSERT_FALSE(ggutils::is_deletion(record1));

    update_record(hdr,2,2400,"CAAAAAAA,CA,C,CAA",record2);
    ASSERT_TRUE(ggutils::is_deletion(record2));
    ASSERT_FALSE(ggutils::is_insertion(record2));
    ASSERT_TRUE( ggutils::bcf1_leq(record1,record2));

    update_record(hdr,2,2400,"G,GGTGTGT",record1);
    update_record(hdr,2,2400,"GGGGTGTGTGT,GGT,G",record2);
    ASSERT_EQ(ggutils::get_variant_rank(record1),1);
    ASSERT_EQ(ggutils::get_variant_rank(record2),2);

    ASSERT_TRUE( ggutils::bcf1_greater_than(record2,record1));
    ASSERT_TRUE( ggutils::bcf1_less_than(record1,record2));
    ASSERT_TRUE(ggutils::is_insertion(record1));
    ggutils::print_variant(record1);
    ggutils::print_variant(record2);

    update_record(hdr,0,11017,"C,T",record1);
    update_record(hdr,1,10000,"T,A",record2);
    ASSERT_FALSE( ggutils::bcf1_equal(record2,record1));
    ASSERT_FALSE( ggutils::bcf1_less_than(record2,record1));
    ASSERT_TRUE( ggutils::bcf1_greater_than(record2,record1));
    bcf_destroy(record1);
    bcf_destroy(record2);
}

// 0/0 0/1 1/1 0/2 1/2 2/2
TEST(UtilTest, GenotypeIndex)
{
    ASSERT_EQ(0, ggutils::get_gl_index(0, 0));
    ASSERT_EQ(1, ggutils::get_gl_index(0, 1));
    ASSERT_EQ(2, ggutils::get_gl_index(1, 1));
    ASSERT_EQ(3, ggutils::get_gl_index(0, 2));
    ASSERT_EQ(4, ggutils::get_gl_index(1, 2));
    ASSERT_EQ(5, ggutils::get_gl_index(2, 2));
}

TEST(UtilTest, phred)
{
    ASSERT_EQ(0, ggutils::phred(1.0));
    ASSERT_EQ(30, ggutils::phred(.001));
    ASSERT_FLOAT_EQ(1.0, ggutils::unphred(0));
    ASSERT_FLOAT_EQ(.001, ggutils::unphred(30));
}


TEST(UtilTest,getploidy)
{
    auto hdr = get_header();
    auto record1 = generate_record(hdr, "chr1\t5420\t.\tC\tA,T,G\t100\tPASS\t.\tGT:GQ:DP:DPF:AD:PL\t1/3:50:16:0:0,12,0,4:396,92,63,368,92,396,276,0,276,285");

    ASSERT_EQ(2,ggutils::get_ploidy(hdr,record1));

    auto record2 = generate_record(hdr, "chr1\t5420\t.\tC\tA,T,G\t100\tPASS\t.\tGT:GQ:DP:DPF:AD\t1:50:16:0:0,12,0,4");
    ASSERT_EQ(1,ggutils::get_ploidy(hdr,record2));
}


TEST(UtilTest,getters)
{
    auto hdr = get_header();
    auto record1 = generate_record(hdr, "chr1\t5420\t.\tC\tA,T,G\t100\tPASS\tMQ=50;AF1000G=.13\tGT:GQ:DP:DPF:AD:PL:SB\t1/3:30:16:0:0,12,0,4:396,92,63,368,92,396,276,0,276,285:1.01");
    int32_t mq,dp;
    ggutils::bcf1_get_one_info_int(hdr,record1,"MQ",mq);
    ggutils::bcf1_get_one_format_int(hdr,record1,"DP",dp);
    ASSERT_EQ(mq,50);
    ASSERT_EQ(dp,16);

    float af;
    ggutils::bcf1_get_one_info_float(hdr,record1,"AF1000G",af);
    ASSERT_FLOAT_EQ(af,.13);
    float sb;
    ggutils::bcf1_get_one_format_float(hdr,record1,"SB",sb);
    ASSERT_FLOAT_EQ(sb,1.01);
}

TEST(UtilTest,bcf1AlleleSwap)
{
    auto hdr = get_header();
    auto record1 = generate_record(hdr, "chr1\t5420\t.\tC\tA,T,G\t100\tPASS\tMQ=50;AF1000G=.13\tGT:GQ:DP:DPF:AD:PL:SB\t1/3:30:16:0:0,12,0,4:396,92,63,368,92,396,276,0,276,285:1.01");

//    ggutils::print_variant(hdr,record1);
    for(int i=1;i<record1->n_allele;i++)
    {
        bcf1_t *record2=bcf_dup(record1);
        ggutils::bcf1_allele_swap(hdr,record2,i,1);
//        ggutils::print_variant(hdr,record2);

        bcf_destroy1(record2);
    }
}