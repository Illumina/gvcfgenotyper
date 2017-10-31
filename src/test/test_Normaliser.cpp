//
// Created by O'Connell, Jared on 10/16/17.
//

#include "gtest/gtest.h"
#include "GVCFReader.hpp"
#include "common.hpp"
#include "test_utils.h"

TEST(Normaliser, mnp_split1)
{
    std::string gvcf_file_name = g_testenv->getBasePath() + "/data/NA12877.tiny.vcf.gz";
    bcf_hdr_t *hdr = bcf_hdr_read(hts_open(gvcf_file_name.c_str(), "r"));
    bcf1_t *record1 = bcf_init1();
    record1->rid = 2;
    record1->pos = 92;
    bcf_update_alleles_str(hdr, record1, "CGG,AGT");
    int32_t ad[2] = {15, 17};
    int32_t pl[3] = {255, 0, 255};
    int32_t gt[2] = {bcf_gt_unphased(0), bcf_gt_unphased(1)};
    int32_t dp = ad[0] + ad[1];
    int32_t dpf = 1;
    float gq = 58;
    bcf_update_format_float(hdr, record1, "GQ", &gq, 1);
    bcf_update_format_int32(hdr, record1, "DPF", &dpf, 1);
    bcf_update_format_int32(hdr, record1, "DP", &dp, 1);
    bcf_update_format_int32(hdr, record1, "AD", &ad, 2);
    bcf_update_format_int32(hdr, record1, "PL", &pl, 3);
    bcf_update_genotypes(hdr, record1, gt, 2);
    vector<bcf1_t *> buffer;
    mnp_split(record1, hdr, buffer);


    htsFile *output_file = hts_open("/dev/null", "wv");
    bcf_hdr_write(output_file, hdr);
    for (auto it = buffer.begin(); it != buffer.end(); it++)
    {
        bcf_write1(output_file, hdr, *it);
    }
    hts_close(output_file);

    bcf1_t *record2 = bcf_init1();
    record2->rid = 2;
    record2->pos = 92;
    bcf_update_alleles_str(hdr, record2, "C,A");
    ASSERT_TRUE(bcf1_equal(record2, buffer[0]));

    bcf1_t *record3 = bcf_init1();
    record3->rid = 2;
    record3->pos = 94;
    bcf_update_alleles_str(hdr, record3, "G,T");
    ASSERT_TRUE(bcf1_equal(record3, buffer[1]));
}

TEST(Normaliser, mnp_split2)
{
    std::string gvcf_file_name = g_testenv->getBasePath() + "/data/NA12877.tiny.vcf.gz";
    bcf_hdr_t *hdr = bcf_hdr_read(hts_open(gvcf_file_name.c_str(), "r"));
    bcf1_t *record1 = bcf_init1();
    record1->rid = 1;
    record1->pos = 100;
    bcf_update_alleles_str(hdr, record1, "CGG,AGG,AGT");
    int32_t ad[3] = {0, 17, 10};
    int32_t pl[6] = {255, 255, 255, 255, 0, 255};
    int32_t gt[2] = {bcf_gt_unphased(1), bcf_gt_unphased(2)};
    int32_t dp = ad[0] + ad[1];
    int32_t dpf = 1;
    float gq = 58;
    bcf_update_format_int32(hdr, record1, "DP", &dp, 1);
    bcf_update_format_float(hdr, record1, "GQ", &gq, 1);
    bcf_update_format_int32(hdr, record1, "DPF", &dpf, 1);
    bcf_update_format_int32(hdr, record1, "AD", &ad, 3);
    bcf_update_format_int32(hdr, record1, "PL", &pl, 6);
    bcf_update_genotypes(hdr, record1, gt, 2);
    vector<bcf1_t *> buffer;
    mnp_split(record1, hdr, buffer);

    htsFile *output_file = hts_open("/dev/null", "wv");
    bcf_hdr_write(output_file, hdr);
    for (auto it = buffer.begin(); it != buffer.end(); it++)
    {
        bcf_write1(output_file, hdr, *it);
    }
    hts_close(output_file);

    bcf1_t *record2 = bcf_init1();
    record2->rid = 1;
    record2->pos = 100;
    bcf_update_alleles_str(hdr, record2, "C,A");
    ASSERT_TRUE(bcf1_equal(record2, buffer[0]));

    bcf1_t *record3 = bcf_init1();
    record3->rid = 1;
    record3->pos = 102;
    bcf_update_alleles_str(hdr, record3, "G,T");
    ASSERT_TRUE(bcf1_equal(record3, buffer[1]));
}

TEST(Normaliser, mnp_split3)
{
    int rid=3;
    int pos=527;
    auto hdr = get_header();
    auto record1 = generate_record(hdr,rid,pos,"CATTCAAGTGC,CGTTTAAGCGA");
    int32_t ad[2] = {10, 17};
    int32_t gt[2] = {bcf_gt_unphased(1), bcf_gt_unphased(2)};
    int32_t dp = ad[0] + ad[1];
    int32_t pl[3] = {255,0,255};
    float gq = 58;
    bcf_update_genotypes(hdr,record1,gt,2);
    bcf_update_format_int32(hdr, record1, "DPI", &dp, 1);
    bcf_update_format_float(hdr, record1, "GQ", &gq, 1);
    bcf_update_format_int32(hdr, record1, "AD", &ad, 2);
    bcf_update_format_int32(hdr, record1, "PL", &pl, 3);
    vector<bcf1_t *> decomposed;
    mnp_split(record1,hdr,decomposed);
    for(auto rec = decomposed.begin();rec!=decomposed.end();rec++)
    {
        Genotype g(hdr,*rec);
        ASSERT_EQ(g.get_dp(),bcf_int32_missing);
    }
}

TEST(Normaliser, unarise1)
{
    std::string ref_file_name = g_testenv->getBasePath() + "/data/tiny.ref.fa";
    std::string gvcf_file_name = g_testenv->getBasePath() + "/data/NA12877.tiny.vcf.gz";
    bcf_hdr_t *hdr = bcf_hdr_read(hts_open(gvcf_file_name.c_str(), "r"));
    Normaliser norm(ref_file_name, hdr);

    bcf1_t *record1 = generate_record(hdr,3,100,"A,C");
    int32_t ad[3] = {11, 19};
    int32_t pl[6] = {255, 0, 255};
    int32_t gt[2] = {bcf_gt_unphased(0), bcf_gt_unphased(1)};
    int32_t dp = ad[0] + ad[1];
    int32_t dpf = 1;
    float gq = 58;
    bcf_update_format_float(hdr, record1, "GQ", &gq, 1);
    bcf_update_format_int32(hdr, record1, "DPF", &dpf, 1);
    bcf_update_format_int32(hdr, record1, "DP", &dp, 1);
    bcf_update_format_int32(hdr, record1, "AD", &ad, 2);
    bcf_update_format_int32(hdr, record1, "PL", &pl, 3);
    bcf_update_genotypes(hdr, record1, gt, 2);
    //print_variant(hdr, record1);


    vector<bcf1_t *> buffer = norm.unarise(record1);
    htsFile *output_file = hts_open("test.out", "wv");
    bcf_hdr_write(output_file, hdr);

    for (auto it = buffer.begin(); it != buffer.end(); it++)
    {
        bcf_write1(output_file, hdr, *it);
        //print_variant(hdr, *it);
    }
    hts_close(output_file);
}


TEST(Normaliser, unarise2)
{
    std::string ref_file_name = g_testenv->getBasePath() + "/data/tiny.ref.fa";
    std::string gvcf_file_name = g_testenv->getBasePath() + "/data/NA12877.tiny.vcf.gz";
    bcf_hdr_t *hdr = bcf_hdr_read(hts_open(gvcf_file_name.c_str(), "r"));
    Normaliser norm(ref_file_name, hdr);

    bcf1_t *record1 = generate_record(hdr,3,100,"A,C,G");
    int32_t ad[3] = {0, 17, 10};
    int32_t pl[6] = {255, 255, 255, 255, 0, 255};
    int32_t gt[2] = {bcf_gt_unphased(1), bcf_gt_unphased(2)};
    int32_t dp = ad[0] + ad[1] + ad[2];
    int32_t dpf = 1;
    float gq = 58;
    bcf_update_format_float(hdr, record1, "GQ", &gq, 1);
    bcf_update_format_int32(hdr, record1, "DPF", &dpf, 1);
    bcf_update_format_int32(hdr, record1, "DP", &dp, 1);
    bcf_update_format_int32(hdr, record1, "AD", &ad, 3);
    bcf_update_format_int32(hdr, record1, "PL", &pl, 6);
    bcf_update_genotypes(hdr, record1, gt, 2);
    //print_variant(hdr, record1);

    vector<bcf1_t *> buffer = norm.unarise(record1);
    htsFile *output_file = hts_open("test.out", "wv");
    bcf_hdr_write(output_file, hdr);

    for (auto it = buffer.begin(); it != buffer.end(); it++)
    {
        bcf_write1(output_file, hdr, *it);
        //print_variant(hdr, *it);
    }
    hts_close(output_file);
}

TEST(Normaliser, unarise3)
{
    std::string ref_file_name = g_testenv->getBasePath() + "/data/tiny.ref.fa";
    std::string gvcf_file_name = g_testenv->getBasePath() + "/data/NA12877.tiny.vcf.gz";
    bcf_hdr_t *hdr = bcf_hdr_read(hts_open(gvcf_file_name.c_str(), "r"));
    Normaliser norm(ref_file_name, hdr);

    bcf1_t *record1 = generate_record(hdr,3,100,"ATT,ATG,CTT");
    int32_t ad[3] = {0, 17, 10};
    int32_t pl[6] = {255, 255, 255, 255, 0, 255};
    int32_t gt[2] = {bcf_gt_unphased(1), bcf_gt_unphased(2)};
    int32_t dp = ad[0] + ad[1];
    int32_t dpf = 1;
    float gq = 58;
    bcf_update_format_float(hdr, record1, "GQ", &gq, 1);
    bcf_update_format_int32(hdr, record1, "DPF", &dpf, 1);
    bcf_update_format_int32(hdr, record1, "DP", &dp, 1);
    bcf_update_format_int32(hdr, record1, "AD", &ad, 3);
    bcf_update_format_int32(hdr, record1, "PL", &pl, 6);
    bcf_update_genotypes(hdr, record1, gt, 2);
    //print_variant(hdr, record1);


    vector<bcf1_t *> buffer = norm.unarise(record1);
    htsFile *output_file = hts_open("test.out", "wv");
    bcf_hdr_write(output_file, hdr);

    for (auto it = buffer.begin(); it != buffer.end(); it++)
    {
        bcf_write1(output_file, hdr, *it);
        //print_variant(hdr, *it);
    }
    hts_close(output_file);
}

//regression test checking that QUAL is correctly propagated by Normaliser::unarise
TEST(Normaliser, unarise4)
{
    int rid=20,pos=83250;
    auto hdr = get_header();
    std::string ref_file_name = g_testenv->getBasePath() + "/data/chr20.100kb.fa";
    Normaliser norm(ref_file_name, hdr);
    auto record1 = generate_record(hdr,rid,pos+1,"TGTGTGTGTG,TC");
    record1->qual = 100;
    int32_t ad[2] = {17, 10};
    int32_t pl[3] = {255, 0,255};
    int32_t gt[2] = {bcf_gt_unphased(0), bcf_gt_unphased(1)};
    float gq = 58;
    bcf_update_format_float(hdr, record1, "GQ", &gq, 1);
    bcf_update_format_int32(hdr, record1, "AD", &ad, 2);
    bcf_update_format_int32(hdr, record1, "PL", &pl, 3);
    bcf_update_genotypes(hdr, record1, gt, 2);
    vector<bcf1_t *> buffer = norm.unarise(record1);
    for (auto it = buffer.begin(); it != buffer.end(); it++)
    {
        print_variant(hdr,*it);
    }
}
//regression test checking that QUAL is correctly propagated by Normaliser::unarise
TEST(Normaliser, unarise5)
{
    int rid = 20, pos = 83250;
    auto hdr = get_header();
    std::string ref_file_name = g_testenv->getBasePath() + "/data/test2/test2.ref.fa";

    Normaliser norm(ref_file_name, hdr);

    auto record1 = generate_record(hdr,"chr1\t5420\t.\tCAAAAAA\tC,A\t423\tPASS\t.\tGT:GQ:AD:PL\t1/2:49:2,15,16:429,123,50,137,0,295");

    vector<bcf1_t *> buffer = norm.unarise(record1);
    auto record2 = generate_record(hdr,"chr1\t5418\t.\tTACAAAA\tT,X\t423\tPASS\t.\tGT\t2/1");
    auto record3 = generate_record(hdr,"chr1	5420	.	CAAAAAA	C,X	423	PASS	.	GT	1/2");

    for (auto it = buffer.begin(); it != buffer.end(); it++)
    {
        print_variant(hdr,*it);
    }
    ASSERT_TRUE(bcf1_equal(record2,buffer[1]));
    ASSERT_TRUE(bcf1_equal(record3,buffer[0]));
}


//regression test checking that QUAL is correctly propagated by Normaliser::unarise
TEST(Normaliser, qual)
{
    int rid=1;
    int pos=4151;
    auto hdr = get_header();
    std::string ref_file_name = g_testenv->getBasePath() + "/data/test2/test2.ref.fa";
    Normaliser norm(ref_file_name, hdr);
    auto record1 = generate_record(hdr,rid,pos+1,"A,ATTT");
    int qual = 221;
    record1->qual = qual;
    int32_t ad[2] = {17, 10};
    int32_t pl[3] = {255, 0,255};
    int32_t gt[2] = {bcf_gt_unphased(0), bcf_gt_unphased(1)};
    float gq = 58;
    bcf_update_format_float(hdr, record1, "GQ", &gq, 1);
    bcf_update_format_int32(hdr, record1, "AD", &ad, 2);
    bcf_update_format_int32(hdr, record1, "PL", &pl, 3);
    bcf_update_genotypes(hdr, record1, gt, 2);
    //print_variant(hdr, record1);
    vector<bcf1_t *> buffer = norm.unarise(record1);
    for (auto it = buffer.begin(); it != buffer.end(); it++)
    {
        ASSERT_FLOAT_EQ((*it)->qual,221);
    }
}
