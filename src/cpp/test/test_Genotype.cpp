//
// Created by O'Connell, Jared on 12/21/17.
//

#include "test_helpers.hh"

extern "C"
{
#include <htslib/vcf.h>
}

#include "GVCFReader.hh"
#include "Genotype.hh"


TEST(Genotype,resolveAlleleConflict)
{
    std::string gvcf_file_name = g_testenv->getBasePath() + "/../test/regression/GG-26.vcf.gz";
    std::string ref_file_name = g_testenv->getBasePath() + "/../test/regression/GG-26.fa";
    int buffer_size = 200;
    GVCFReader reader(gvcf_file_name, ref_file_name, buffer_size);
    bcf_hdr_t *hdr = reader.GetHeader();
    multiAllele m;
    m.Init(hdr);
    m.SetPosition(reader.Front()->rid, reader.Front()->pos);
    auto variants = reader.GetAllVariantsInInterval(reader.Front()->rid, reader.Front()->pos);
    ASSERT_TRUE(variants.first!=variants.second);
    for (auto rec = variants.first; rec != variants.second; rec++)
        m.Allele(*rec);
    auto sample_variants = reader.GetAllVariantsUpTo(m.GetMax());
    Genotype g(hdr,sample_variants,m);
    for(int i=0;i< g.num_allele();i++)
    {
//        std::cerr << g.get_ad(i) << " " << g.get_adf(i)<< "+"<<g.get_adr(i)<<"="<<g.get_adf(i)+g.get_adf(i)<<std::endl;
        ASSERT_EQ(g.ad(i), g.adf(i)+ g.adr(i));
    }
}

//regression test checking that Genotype correctly propagates FORMAT fields
TEST(Genotype,format)
{
    int rid=1;
    int pos=97473;
    auto hdr = get_header();
    std::string ref_file_name = g_testenv->getBasePath() + "/../test/test2/test2.ref.fa";
    Normaliser norm(ref_file_name, hdr);
    auto record1 = generate_record(hdr,rid,pos+1,"G,GA,GAA");
    int qual = 786;
    record1->qual = qual;
    int32_t ad[3] = {2,24,17};
    int32_t adf[3] = {2,12,10};
    int32_t adr[3] = {0,12,7};
    int32_t pl[6] = {405,132,56,188,0,121};
    int32_t gt[2] = {bcf_gt_unphased(1), bcf_gt_unphased(2)};
    float gq = 58;
    bcf_update_format_float(hdr, record1, "GQ", &gq, 1);
    bcf_update_format_int32(hdr, record1, "AD", &ad, 3);
    bcf_update_format_int32(hdr, record1, "ADR", &adr, 3);
    bcf_update_format_int32(hdr, record1, "ADF", &adf, 3);
    bcf_update_format_int32(hdr, record1, "PL", &pl, 6);
    bcf_update_genotypes(hdr, record1, gt, 2);
//    print_variant(hdr, record1);
    vector<bcf1_t *> buffer;
    norm.Unarise(record1, buffer);
    size_t idx = 0;
    for (auto it = buffer.begin(); it != buffer.end(); it++)
    {
        std::cerr<<"idx="<<idx<<std::endl;
        ggutils::print_variant(hdr,*it);
        Genotype g(hdr,*it);
        ASSERT_FLOAT_EQ(g.gq(),gq);

        ASSERT_EQ(g.ad(0),ad[0]);
        ASSERT_EQ(g.ad(1),ad[1+idx]);

        ASSERT_EQ(g.adf(0),adf[0]);
        ASSERT_EQ(g.adf(1),adf[1+idx]);

        ASSERT_EQ(g.adr(0),adr[0]);
        ASSERT_EQ(g.adr(1),adr[1+idx]);

        ++idx;
    }
}

TEST(Genotype,CallGenotypeDiploid)
{
auto hdr = get_header();
    auto record1 = generate_record(hdr, "chr1\t10002558\t.\tG\tA\t499\tPASS\tMQ=60\tGT:GQ:GQX:DP:DPF:AD:ADF:ADR:SB:FT:PL\t0/0:102:30:35:0:0,35:0,15:0,20:-55.2:PASS:370,105,0");
    Genotype g1(hdr,record1);
    ASSERT_EQ(bcf_gt_allele(g1.gt(0)),0);
    ASSERT_EQ(bcf_gt_allele(g1.gt(1)),0);
    g1.CallGenotype();
    ASSERT_EQ(bcf_gt_allele(g1.gt(0)),1);
    ASSERT_EQ(bcf_gt_allele(g1.gt(1)),1);

    auto record2 = generate_record(hdr, "chr1\t10002558\t.\tG\tA\t499\tPASS\tMQ=60\tGT:GQ:GQX:DP:DPF:AD:ADF:ADR:SB:FT:PL\t0/0:102:30:35:0:0,35:0,15:0,20:-55.2:PASS:370,0,105");
    Genotype g2(hdr,record2);
    ASSERT_EQ(bcf_gt_allele(g2.gt(0)),0);
    ASSERT_EQ(bcf_gt_allele(g2.gt(1)),0);
    g2.CallGenotype();
    ASSERT_EQ(bcf_gt_allele(g2.gt(0)),0);
    ASSERT_EQ(bcf_gt_allele(g2.gt(1)),1);

    auto record3 = generate_record(hdr, "chr1\t10002558\t.\tG\tA\t499\tPASS\tMQ=60\tGT:GQ:GQX:DP:DPF:AD:ADF:ADR:SB:FT:PL\t1/1:102:30:35:0:0,35:0,15:0,20:-55.2:PASS:0,370,105");
    Genotype g3(hdr,record3);
    ASSERT_EQ(bcf_gt_allele(g3.gt(0)),1);
    ASSERT_EQ(bcf_gt_allele(g3.gt(1)),1);
    g3.CallGenotype();
    ASSERT_EQ(bcf_gt_allele(g3.gt(0)),0);
    ASSERT_EQ(bcf_gt_allele(g3.gt(1)),0);

    auto record4 = generate_record(hdr,"chr1\t10587268\t.\tC\tCAAA,CAA\t276\tPASS\tMQ=58\tGT:GQ:GQX:DPI:AD:ADF:ADR:FT:PL\t0/0:28:11:27:1,15,9:0,9,4:1,6,5:PASS:288,79,31,144,0,103");
    Genotype g4(hdr,record4);
    ASSERT_EQ(bcf_gt_allele(g4.gt(0)),0);
    ASSERT_EQ(bcf_gt_allele(g4.gt(1)),0);
    g4.CallGenotype();
    ASSERT_EQ(bcf_gt_allele(g4.gt(0)),1);
    ASSERT_EQ(bcf_gt_allele(g4.gt(1)),2);
}

TEST(Genotype,CallGenotypeHaploid)
{
    auto hdr = get_header();
    auto record1 = generate_record(hdr, "chrX\t10027690\t.\tG\tA\t315\tPASS\tMQ=60\tGT:GQ:GQX:DP:DPF:AD:ADF:ADR:SB:FT:PL\t0:342:21:20:2:0,20:0,13:0,7:-30.7:PASS:349,0");
    Genotype g1(hdr, record1);
    ASSERT_EQ(bcf_gt_allele(g1.gt(0)), 0);
    g1.CallGenotype();
    ASSERT_EQ(bcf_gt_allele(g1.gt(0)), 1);

    auto record2 = generate_record(hdr, "chrX\t52891584\t.\tC\tA,T\t282\tPASS\tMQ=60\tGT:GQ:GQX:DP:DPF:AD:ADF:ADR:SB:FT:PL\t1:35:30:32:1:0,16,16:0,5,10:0,11,6:-30.6:PASS:316,35,0");
    Genotype g2(hdr, record2);
    ASSERT_EQ(bcf_gt_allele(g2.gt(0)), 1);
    g2.CallGenotype();
    ASSERT_EQ(bcf_gt_allele(g2.gt(0)), 2);
}
