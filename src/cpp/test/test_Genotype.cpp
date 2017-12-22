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
    std::string gvcf_file_name = g_testenv->getBasePath() + "/../test/regression/GG-19.vcf.gz";
    std::string ref_file_name = g_testenv->getBasePath() + "/../test/regression/GG-19.fa";
    int buffer_size = 200;
    GVCFReader reader(gvcf_file_name, ref_file_name, buffer_size);
    bcf_hdr_t *hdr = reader.get_header();
    multiAllele m;
    m.init(hdr);
    m.setPosition(reader.front()->rid,reader.front()->pos);
    auto variants = reader.get_all_variants_in_interval(reader.front()->rid,reader.front()->pos);
    ASSERT_TRUE(variants.first!=variants.second);
    for (auto rec = variants.first; rec != variants.second; rec++)
        m.allele(*rec);
    auto sample_variants = reader.get_all_variants_up_to(m.get_max());
    Genotype g(hdr,sample_variants,m);
    for(int i=0;i<g.get_num_allele();i++)
    {
//        std::cerr << g.get_ad(i) << " " << g.get_adf(i)<< "+"<<g.get_adr(i)<<"="<<g.get_adf(i)+g.get_adf(i)<<std::endl;
        ASSERT_EQ(g.get_ad(i),g.get_adf(i)+g.get_adr(i));
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
    norm.unarise(record1,buffer);
    size_t idx = 0;
    for (auto it = buffer.begin(); it != buffer.end(); it++)
    {
        std::cerr<<"idx="<<idx<<std::endl;
        ggutils::print_variant(hdr,*it);
        Genotype g(hdr,*it);
        ASSERT_FLOAT_EQ(g.get_gq(),gq);

        ASSERT_EQ(g.get_ad(0),ad[0]);
        ASSERT_EQ(g.get_ad(1),ad[1+idx]);

        ASSERT_EQ(g.get_adf(0),adf[0]);
        ASSERT_EQ(g.get_adf(1),adf[1+idx]);

        ASSERT_EQ(g.get_adr(0),adr[0]);
        ASSERT_EQ(g.get_adr(1),adr[1+idx]);

        ++idx;
    }
}
