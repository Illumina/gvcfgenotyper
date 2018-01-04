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
    ASSERT_EQ(db.dp(), 20);
    buf.interpolate(0, 99, 99, db);
    ASSERT_EQ(db.dp(), 20);
    buf.interpolate(0, 100, 100, db);
    ASSERT_EQ(db.dp(), 30);
    buf.interpolate(0, 95, 104, db);
    ASSERT_EQ(db.dp(), 25);
    buf.interpolate(0, 95, 114, db);
    ASSERT_EQ(db.dp(), 30);
    buf.interpolate(0, 95, 154, db);
    ASSERT_EQ(db.dp(), 37);
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

//This tests a fairly tricky case of two different insertions starting at the same position within the same sample.
TEST(VariantBuffer, test2)
{
    multiAllele m;
    auto hdr = get_header();
    m.init(hdr);
    std::string ref_file_name = g_testenv->getBasePath() + "/../test/test2/test2.ref.fa";
    Normaliser norm(ref_file_name, hdr);
    auto record1 = generate_record(hdr,"chr1\t7832\trs112070696\tC\tCTAAATAAATAAA,CTAAATAAATAAATAAA\t559\tPASS\t.\tGT:GQ:GQX:DPI:AD:ADF:ADR:FT:PL\t1/2:150:15:42:0,11,11:0,4,5:0,7,6:PASS:601,226,169,225,0,169");
    m.setPosition(record1->rid,record1->pos);
    vector<bcf1_t *> buffer;
    norm.unarise(record1,buffer);
    VariantBuffer v;
    int count=0;
    for (auto it = buffer.begin(); it != buffer.end(); it++)
    {
        v.push_back(hdr,*it);
        ASSERT_EQ(m.allele(*it),++count);
    }
    v.flush_buffer(m.get_max());
    ASSERT_TRUE(v.empty());
}

TEST(GVCFReader, readMNP)
{
    std::string gvcf_file_name = g_testenv->getBasePath() + "/../test/mnp.genome.vcf";
    std::string ref_file_name = g_testenv->getBasePath() + "/../test/tiny.ref.fa";
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
    std::string expected_output_file = g_testenv->getBasePath() + "/../test/NA12877.tiny.vcf.gz.expected";
    std::string gvcf_file_name = g_testenv->getBasePath() + "/../test/NA12877.tiny.vcf.gz";
    std::string ref_file_name = g_testenv->getBasePath() + "/../test/tiny.ref.fa";
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
//        print_variant(hdr,line);
        ofs << bcf_hdr_id2name(hdr, line->rid) << ":" << line->pos + 1 << ":" << line->d.allele[0] << ":" << line->d.allele[1] << std::endl;
        if (ggutils::is_snp(line))
        {
            reader.get_depth(line->rid, line->pos, line->pos, db);
            if (bcf_get_format_int32(hdr, line, "DP", &dp, &nval) == 1)
            {
                ASSERT_EQ(db.dp(), *dp);
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


