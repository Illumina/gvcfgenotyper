#include "test_helpers.hh"

#include "GVCFMerger.hh"
#include "StringUtil.hh"


TEST(multiAllele,test1)
{
    int rid=1;
    int pos=99;
    auto hdr = get_header();
    auto rec1 = generate_record(hdr,rid,pos,"C,G");
    auto rec2 = generate_record(hdr,rid,pos,"C,A");
    auto rec3 = generate_record(hdr,rid,200,"CTG,C");
    auto rec4 = generate_record(hdr,rid,pos,"CTGG,C,CAAAAAAAATGG");
    auto rec5 = generate_record(hdr,rid,pos,"CTGG,CAAAAAAAATGG,C");

    multiAllele m;
    m.init(hdr);
    m.setPosition(rid,pos-1);
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
    ggutils::print_variant(hdr,v);

    m.get_max();

    auto truth = generate_record(hdr,rid,pos,"CTGG,GTGG,ATGG,C,CAAAAAAAATGG");     //chr1:100:CTGG:GTGG,ATGG,C,CAAAAAAAATGG
    ASSERT_TRUE( ggutils::bcf1_all_equal(truth,v));
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
    std::string test_base = g_testenv->getBasePath() + "/../test/test2/";
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

    std::string ref_file_name = g_testenv->getBasePath() + "/../test/test2/test2.ref.fa";
    std::string output_file_name = "test.out";//std::tmpnam(NULL);
    std::cerr << "Outputting to " << output_file_name << std::endl;
    GVCFMerger g(files, output_file_name, "z", ref_file_name, buffer_size);
    g.write_vcf();
}

TEST(GVCFMerger, likelihood)
{
    auto hdr = get_header();
    std::string ref_file_name = g_testenv->getBasePath() + "/../test/test2/test2.ref.fa";
    Normaliser norm(ref_file_name, hdr);
    auto record1 = generate_record(hdr,"chr1\t85677\t.\tT\tA\t191\tPASS\t.\tGT:GQ:GQX:DP:DPF:AD:ADF:ADR:SB:FT:PL\t0/1:224:30:42:0:22,20:13,11:9,9:-26:PASS:226,0,257");

    std::cerr <<"Input:"<<std::endl;
    ggutils::print_variant(hdr,record1);
    multiAllele m;
    m.init(hdr);
    m.setPosition(record1->rid,record1->pos);
    m.allele(record1);

    ggutils::print_variant(hdr,record1);

    ggutils::vcf_data_t d(2,2,2);
    Genotype g(hdr,record1);
    g.propagate_format_fields(0,2,&d);
    ASSERT_EQ(d.pl[0],226);
    ASSERT_EQ(d.pl[1],0);
    ASSERT_EQ(d.pl[2],257);

    g.propagate_format_fields(1,2,&d);
    ASSERT_EQ(d.pl[3],226);
    ASSERT_EQ(d.pl[4],0);
    ASSERT_EQ(d.pl[5],257);

}
