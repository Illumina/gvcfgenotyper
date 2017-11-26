#include "test_helpers.hh"

extern "C"
{
#include <htslib/vcf.h>
}

#include "GVCFReader.hh"



TEST(DepthBlock, intersects)
{
    DepthBlock db1(0, 100, 199, 20, 1, 30,2);
    DepthBlock db2(0, 100, 100, 20, 1, 30,2);
    ASSERT_EQ(db1.intersect_size(db2), 1);
    ASSERT_EQ(db1.intersect_size(db2), db2.intersect_size(db1));

    DepthBlock db3(0, 300, 401, 20, 1, 30,2);
    ASSERT_EQ(db1.intersect_size(db3), 0);

    DepthBlock db4(1, 300, 401, 20, 1, 30,2);
    ASSERT_EQ(db1.intersect_size(db4), 0);

    DepthBlock db5(0, 190, 401, 20, 1, 30,2);
    ASSERT_EQ(db1.intersect_size(db5), 10);

    DepthBlock db6(0, 100, 100, 20, 1, 30,2);
    DepthBlock db7(0, 100, 100, 20, 1, 30,2);
    ASSERT_EQ(db6.intersect_size(db7), 1);

    DepthBlock db8(1, 10000, 10000, 20, 1, 30,2);
    DepthBlock db9(1, 10000, 10000, 20, 1, 30,2);
    ASSERT_EQ(db8.intersect_size(db9), 1);

}

TEST(DepthBuffer, interpolate)
{
    DepthBuffer buf;
    buf.push_back(DepthBlock(0, 0, 99, 20, 1, 30,2));
    buf.push_back(DepthBlock(0, 100, 109, 30, 1, 30,2));
    buf.push_back(DepthBlock(0, 110, 200, 40, 1, 30,2));
    DepthBlock db;
    buf.interpolate(0, 90, 95, db);
    ASSERT_EQ(db._dp, 20);
    buf.interpolate(0, 99, 99, db);
    ASSERT_EQ(db._dp, 20);
    buf.interpolate(0, 100, 100, db);
    ASSERT_EQ(db._dp, 30);
    buf.interpolate(0, 95, 104, db);
    ASSERT_EQ(db._dp, 25);
    buf.interpolate(0, 95, 114, db);
    ASSERT_EQ(db._dp, 30);
    buf.interpolate(0, 95, 154, db);
    ASSERT_EQ(db._dp, 37);
}

TEST(VariantBuffer, test1)
{
    int rid=1;
    int pos=99;
    auto hdr = get_header();
    auto rec1 = generate_record(hdr,rid,pos,"C,G");

    VariantBuffer v;
    v.push_back(hdr,rec1);
    v.flush_buffer(rec1);
    ASSERT_TRUE(v.empty());

    auto rec2 = generate_record(hdr,rid,pos,"C,G");
    auto rec3 = generate_record(hdr,rid,pos,"C,CG");
    v.push_back(hdr,bcf_dup(rec2));
    v.push_back(hdr,bcf_dup(rec3));

    v.flush_buffer(rec2);
    ASSERT_FALSE(v.empty());

    v.flush_buffer(rec3);
    ASSERT_TRUE(v.empty());
}

TEST(GVCFReader, readMNP)
{
    std::string gvcf_file_name = g_testenv->getBasePath() + "/data/mnp.genome.vcf";
    std::string ref_file_name = g_testenv->getBasePath() + "/data/tiny.ref.fa";
    GVCFReader reader(gvcf_file_name, ref_file_name, 1000);
    const bcf_hdr_t *hdr = reader.get_header();
    bcf1_t *line = reader.pop();
    int32_t *dp = nullptr, nval = 0;
    while (line != nullptr)
    {
        if (ggutils::is_snp(line))
        {
            ASSERT_EQ(bcf_get_format_int32(hdr, line, "DP", &dp, &nval),-3);
        }
        bcf_destroy(line);
        line = reader.pop();
    }
    free(dp);
}

// //GVCFReader should match this VID output
// //bcftools norm -m -any data/NA12877.tiny.vcf.gz | bcftools norm -f data/tiny.ref.fa | bcftools query -i 'ALT!="."' -f '%CHROM:%POS:%REF:%ALT\n' > data/NA12877.tiny.vcf.gz.expected 
TEST(GVCFReader, readAGVCF)
{
    std::string expected_output_file = g_testenv->getBasePath() + "/data/NA12877.tiny.vcf.gz.expected";
    std::string gvcf_file_name = g_testenv->getBasePath() + "/data/NA12877.tiny.vcf.gz";
    std::string ref_file_name = g_testenv->getBasePath() + "/data/tiny.ref.fa";
    char tn[] = "/tmp/tmpvcf-XXXXXXX";
    int fd = mkstemp(tn);
    if (fd < 1)
    {
        error("Cannot create temp file: %s", strerror(errno));
    }
    close(fd);

    std::ofstream ofs(tn, std::ofstream::out);
    int buffer_size = 200;
    GVCFReader reader(gvcf_file_name, ref_file_name, buffer_size);
    const bcf_hdr_t *hdr = reader.get_header();
    bcf1_t *line = reader.pop();
    int32_t *dp = nullptr, nval = 0;
    DepthBlock db;
    while (line != nullptr)
    {
//        ggutils::print_variant(hdr,line);
        ofs << bcf_hdr_id2name(hdr, line->rid) << ":" << line->pos + 1 << ":" << line->d.allele[0] << ":" << line->d.allele[1] << std::endl;
        if (ggutils::is_snp(line))
        {
            reader.get_depth(line->rid, line->pos, line->pos, db);
            if (bcf_get_format_int32(hdr, line, "DP", &dp, &nval) == 1)
            {
                ASSERT_EQ(db._dp, *dp);
            }
        }
        if(line->n_allele>1)
        {
            float gq;
            int ngq=0;
            if(bcf_get_format_int32(reader.get_header(),line,"GQ",&gq,&ngq)==1)
            {
                Genotype g(reader.get_header(),line);
                ASSERT_EQ(g.get_gq(),(int)gq);
            }
        }
        bcf_destroy(line);
        line = reader.pop();
    }
    free(dp);
    ofs.close();
    ASSERT_TRUE(reader.empty());
    const std::string diffcmd = std::string("diff -I '^#' ") + tn + " " + expected_output_file;
    int r = system(diffcmd.c_str());
    if (r != 0)
    {
        error("Difference detected in test case %s: %s\n\n", gvcf_file_name.c_str(), diffcmd.c_str());
    }
    else
    {
        remove(tn);
    }
}

//regression test checking that Genotype correctly propagates FORMAT fields
TEST(Genotype,format)
{
    int rid=1;
    int pos=97473;
    auto hdr = get_header();
    std::string ref_file_name = g_testenv->getBasePath() + "/data/test2/test2.ref.fa";
    Normaliser norm(ref_file_name, hdr);
    auto record1 = generate_record(hdr,rid,pos+1,"G,GA,GAA");
    int qual = 786;
    record1->qual = qual;
    int32_t ad[3] = {2,24,17};
    int32_t pl[6] = {405,132,56,188,0,121};
    int32_t gt[2] = {bcf_gt_unphased(1), bcf_gt_unphased(2)};
    float gq = 58;
    bcf_update_format_float(hdr, record1, "GQ", &gq, 1);
    bcf_update_format_int32(hdr, record1, "AD", &ad, 3);
    bcf_update_format_int32(hdr, record1, "PL", &pl, 6);
    bcf_update_genotypes(hdr, record1, gt, 2);
//    ggutils::print_variant(hdr, record1);
    vector<bcf1_t *> buffer;
    norm.unarise(record1,buffer);
    for (auto it = buffer.begin(); it != buffer.end(); it++)
    {
//        ggutils::print_variant(hdr,*it);
        Genotype g(hdr,*it);
        ASSERT_FLOAT_EQ(g.get_gq(),gq);
    }
}

//tests marginalisation and propagation of likelihoods in our genotype class.
TEST(Genotype,likelihood1)
{
    auto hdr = get_header();
    auto record1 = generate_record(hdr, "chr1\t5420\t.\tC\tA,T\t100\tPASS\t.\tGT:GQ:DP:DPF:AD:PL\t1/2:50:25:0:1,13,11:526,230,276,214,0,266");
    Genotype g(hdr, record1);
    g.print();
    Genotype g1 = g.marginalise(1);
    g1.print();
    Genotype g2 = g.marginalise(2);
    g2.print();
    ASSERT_EQ(0,g1.get_pl(1,2));
    ASSERT_EQ(0,g2.get_pl(1,2));
}

TEST(Genotype,likelihood2)
{
    auto hdr = get_header();
    auto record1 = generate_record(hdr,"chr1\t5420\t.\tC\tA,T,G\t100\tPASS\t.\tGT:GQ:DP:DPF:AD:PL\t1/3:50:16:0:0,12,0,4:396,92,63,368,92,396,276,0,276,285");
    Genotype g(hdr,record1);
    g.print();
    Genotype g1 = g.marginalise(1);
    g1.print();
    Genotype g2 = g.marginalise(2);
    g2.print();
    Genotype g3 = g.marginalise(3);
    g3.print();
    ASSERT_EQ(0,g1.get_pl(1,2));
    ASSERT_EQ(0,g2.get_pl(2,2));
    ASSERT_EQ(0,g3.get_pl(1,2));
}
