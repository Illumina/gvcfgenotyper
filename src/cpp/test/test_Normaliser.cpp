//
// Created by O'Connell, Jared on 10/16/17.
//

#include "test_helpers.hh"
#include "GVCFReader.hh"

//regression test checking that QUAL is correctly propagated by Normaliser::unarise
TEST(Normaliser, qual)
{
    int rid=1;
    int pos=4151;
    auto hdr = get_header();
    std::string ref_file_name = g_testenv->getBasePath() + "/../test/test2/test2.ref.fa";
    Normaliser norm(ref_file_name);
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
    vector<bcf1_t *> buffer;
    norm.Unarise(record1, buffer,hdr);
    for (auto it = buffer.begin(); it != buffer.end(); it++)
    {
        ASSERT_FLOAT_EQ((*it)->qual,221);
    }
}

TEST(Normaliser, mnp_decompose1)
{
    std::string gvcf_file_name = g_testenv->getBasePath() + "/../test/NA12877.tiny.vcf.gz";
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
    mnp_decompose(record1, hdr, buffer);


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
    ASSERT_TRUE( ggutils::bcf1_equal(record2, buffer[0]));

    bcf1_t *record3 = bcf_init1();
    record3->rid = 2;
    record3->pos = 94;
    bcf_update_alleles_str(hdr, record3, "G,T");
    ASSERT_TRUE( ggutils::bcf1_equal(record3, buffer[1]));
}

TEST(Normaliser, mnp_decompose2)
{
    std::string gvcf_file_name = g_testenv->getBasePath() + "/../test/NA12877.tiny.vcf.gz";
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
    mnp_decompose(record1, hdr, buffer);

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
    ASSERT_TRUE(ggutils::bcf1_equal(record2, buffer[0]));

    bcf1_t *record3 = bcf_init1();
    record3->rid = 1;
    record3->pos = 102;
    bcf_update_alleles_str(hdr, record3, "G,T");
    ASSERT_TRUE( ggutils::bcf1_equal(record3, buffer[1]));
}

TEST(Normaliser, mnp_decompose3)
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
    mnp_decompose(record1,hdr,decomposed);
    for(auto rec = decomposed.begin();rec!=decomposed.end();rec++)
    {
        Genotype g(hdr,*rec);
//        ASSERT_EQ(g.get_dp(),bcf_int32_missing);
    }
}

TEST(Normaliser, unarise1)
{
    std::string ref_file_name = g_testenv->getBasePath() + "/../test/tiny.ref.fa";
    std::string gvcf_file_name = g_testenv->getBasePath() + "/../test/NA12877.tiny.vcf.gz";
    bcf_hdr_t *hdr = bcf_hdr_read(hts_open(gvcf_file_name.c_str(), "r"));
    Normaliser norm(ref_file_name);

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


    vector<bcf1_t *> buffer;
    norm.Unarise(record1, buffer,hdr);
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
    std::string ref_file_name = g_testenv->getBasePath() + "/../test/tiny.ref.fa";
    std::string gvcf_file_name = g_testenv->getBasePath() + "/../test/NA12877.tiny.vcf.gz";
    bcf_hdr_t *hdr = bcf_hdr_read(hts_open(gvcf_file_name.c_str(), "r"));
    Normaliser norm(ref_file_name);

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

    vector<bcf1_t *> buffer;
    norm.Unarise(record1, buffer,hdr);
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
    std::string ref_file_name = g_testenv->getBasePath() + "/../test/tiny.ref.fa";
    std::string gvcf_file_name = g_testenv->getBasePath() + "/../test/NA12877.tiny.vcf.gz";
    bcf_hdr_t *hdr = bcf_hdr_read(hts_open(gvcf_file_name.c_str(), "r"));
    Normaliser norm(ref_file_name);

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


    vector<bcf1_t *> buffer;
    norm.Unarise(record1, buffer,hdr);
    htsFile *output_file = hts_open("test.out", "wv");
    bcf_hdr_write(output_file, hdr);

    for (auto it = buffer.begin(); it != buffer.end(); it++)
    {
        bcf_write1(output_file, hdr, *it);
        //print_variant(hdr, *it);
    }
    hts_close(output_file);
}


TEST(Normaliser, unarise4)
{
    int rid=20,pos=83250;
    auto hdr = get_header();
    std::string ref_file_name = g_testenv->getBasePath() + "/../test/chr20.100kb.fa";
    Normaliser norm(ref_file_name);
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
    vector<bcf1_t *> buffer;
    norm.Unarise(record1, buffer,hdr);
    for (auto it = buffer.begin(); it != buffer.end(); it++)
    {
        ggutils::print_variant(hdr,*it);
    }
}


TEST(Normaliser, unarise6)
{
    auto hdr = get_header();
    std::string ref_file_name = g_testenv->getBasePath() + "/../test/test2/test2.ref.fa";
    Normaliser norm(ref_file_name);
    auto record1 = generate_record(hdr,"chr1\t5420\t.\tCAAAAAA\tC,A\t423\tPASS\t.\tGT:GQ:GQX:DPI:AD:ADF:ADR:FT:PL\t1/2:49:7:54:2,15,16:1,12,0:1,3,16:PASS:429,123,50,137,0,295");
    multiAllele m;
    m.Init(hdr);
    m.SetPosition(record1->rid, record1->pos);

    std::cerr <<"Input:"<<std::endl;
    ggutils::print_variant(hdr,record1);
    vector<bcf1_t *> buffer;
    norm.Unarise(record1, buffer,hdr);
    std::cerr <<"Output:"<<std::endl;
    for (auto it = buffer.begin(); it != buffer.end(); it++)
    {
        ggutils::print_variant(hdr,*it);
        m.Allele(*it);
    }
//    ggutils::print_variant(hdr,m.get_max());
}

TEST(Normaliser, unarise7)
{
    auto hdr = get_header();
    std::string ref_file_name = g_testenv->getBasePath() + "/../test/test2/test2.ref.fa";
    Normaliser norm(ref_file_name);
    auto record1 = generate_record(hdr,"chr1\t95593\t.\tAAAAAAG\tAAA,A\t738\tLowGQX\t.\tGT:GQ:GQX:DPI:AD:ADF:ADR:FT:PL\t1/2:151:0:29:1,15,14:1,7,8:0,8,6:LowGQX:816,279,187,308,0,231");
    std::cerr <<"Input:"<<std::endl;
    ggutils::print_variant(hdr,record1);
    vector<bcf1_t *> buffer;
    norm.Unarise(record1, buffer,hdr);
    std::cerr <<"Output:"<<std::endl;
    for (auto it = buffer.begin(); it != buffer.end(); it++)
    {
        ggutils::print_variant(hdr,*it);
    }
}

TEST(Normaliser, unarise8)
{
    auto hdr = get_header();
    std::string ref_file_name = g_testenv->getBasePath() + "/../test/test2/test2.ref.fa";
    Normaliser norm(ref_file_name);
    auto record1 = generate_record(hdr,"chr1\t62430033\t.\tC\tT,A\t391\tPASS\t.\tGT:GQ:GQX:DP:DPF:AD:ADF:ADR:SB:FT:PL\t1/2:142:30:28:2:0,16,12:0,7,3:0,9,9:-41.5:PASS:370,214,166,242,0,206");
    std::cerr <<"Input:"<<std::endl;
    ggutils::print_variant(hdr,record1);
    vector<bcf1_t *> buffer;
    norm.Unarise(record1, buffer,hdr);
    std::cerr <<"Output:"<<std::endl;
    for (auto it = buffer.begin(); it != buffer.end(); it++)
    {
        ggutils::print_variant(hdr,*it);
    }
}

TEST(Normaliser, unarise9)
{
    auto hdr = get_header();
    std::string ref_file_name = g_testenv->getBasePath() + "/../test/test2/test2.ref.fa";
    Normaliser norm(ref_file_name);
    auto record1 = generate_record(hdr,"chr1\t1\t.\tA\tT,G\t0\tLowGQX\tSNVHPOL=2;MQ=33\tGT:GQ:GQX:DP:DPF:AD:ADF:ADR:SB:FT:PL\t0/2:9:0:9:0:7,1,1:5,1,1:2,0,0:0:LowGQX:14,11,149,0,117,138");
    multiAllele m;
    m.Init(hdr);
    m.SetPosition(record1->rid, record1->pos);

    std::cerr <<"Input:"<<std::endl;
    ggutils::print_variant(hdr,record1);
    std::vector<bcf1_t *> buffer;
    norm.Unarise(record1, buffer,hdr);
    std::cerr <<"Output:"<<std::endl;
    std::deque<bcf1_t *> q;
    for (auto it = buffer.begin(); it != buffer.end(); it++)
    {
        ggutils::print_variant(hdr,*it);
        m.Allele(*it);
        q.push_back(*it);
    }
    std::cerr<<std::endl;
    pair<std::deque<bcf1_t *>::iterator,std::deque<bcf1_t *>::iterator> i(q.begin(),q.end());
    Genotype g(hdr,i,m);
}

TEST(Normaliser, unarise10)
{
    auto hdr = get_header();
    std::string ref_file_name = g_testenv->getBasePath() + "/../test/test2/test2.ref.fa";
    Normaliser norm(ref_file_name);

    auto record1 = generate_record(hdr,"chr1\t4\trs863224454;rs863224781\tGCA\tCTG,CCA\t2\tLowGQX\tRU=.,.;REFREP=.,.;IDREP=.,.;MQ=60;OLD_VARIANT=chr5:148407707:GACG/GCTG/GCCA;clinvar=1|likely_pathogenic,2|uncertain_significance;GMAF=A|0.009385,A|0.009385;CSQT=1|SH3TC2|NM_024577.3|missense_variant,2|SH3TC2|NM_024577.3|missense_variant\tGT:GQ:GQX:DPI:AD:ADF:ADR:FT:PL\t0/1:5:0:56:4,3,3:1,1,1:3,2,2:LowGQX:43,8,8,16,0,29");
    std::vector<bcf1_t *> buffer;
    norm.Unarise(record1, buffer,hdr);
    std::cerr << "Output:" << std::endl;
    for (auto it = buffer.begin(); it != buffer.end(); it++)
    {
        Genotype g(hdr,*it);
        ASSERT_TRUE(g.ad(0)>=0);
        ASSERT_TRUE(g.ad(1)>=0);
        ggutils::print_variant(hdr, *it);
    }
}

TEST(Normaliser, CollapseRecords1)
{
    auto hdr = get_header();
    std::string ref_file_name = g_testenv->getBasePath() + "/../test/test2/test2.ref.fa";
    Normaliser norm(ref_file_name);
    auto record1 = generate_record(hdr,"chr1\t1\t.\tA\tT,G\t0\tLowGQX\tSNVHPOL=2;MQ=33\tGT:GQ:GQX:DP:DPF:AD:ADF:ADR:SB:FT:PL\t0/2:9:0:9:0:7,1,1:5,1,1:2,0,0:0:LowGQX:14,11,149,0,117,138");
    multiAllele m;
    m.Init(hdr);
    m.SetPosition(record1->rid, record1->pos);

    //std::cerr <<"Input:"<<std::endl;    ggutils::print_variant(hdr,record1);
    std::vector<bcf1_t *> buffer;
    norm.Unarise(record1, buffer,hdr);
//    std::cerr <<"Output:"<<std::endl;
    std::deque<bcf1_t *> q;
    for (auto it = buffer.begin(); it != buffer.end(); it++)
    {
  //      ggutils::print_variant(hdr,*it);
        m.Allele(*it);
        q.push_back(*it);
    }
    pair<std::deque<bcf1_t *>::iterator,std::deque<bcf1_t *>::iterator> i(q.begin(),q.end());
    auto collapsed_record = CollapseRecords(hdr,i);
//    ggutils::print_variant(hdr,collapsed_record);
    ASSERT_TRUE(ggutils::bcf1_all_equal(collapsed_record,record1));
}

TEST(Normaliser, CollapseRecords2)
{
    auto hdr = get_header();
    auto rec1 = generate_record(hdr,"chr1\t6700131\t.\tC\tCA\t186\tPASS\tMQ=60\tGT:GQ:GQX:DPI:AD:ADF:ADR:FT:PL\t0/1:191:13:41:20,18:10,8:10,10:PASS:189,0,197");
    auto rec2 = generate_record(hdr,"chr1\t6700131\t.\tC\tCAA\t0\tPASS\tMQ=60\tGT:GQ:GQX:DPI:AD:ADF:ADR:FT:PL\t0/0:105:105:41:37,0:18,0:19,0:PASS:0,108,561");
    auto rec3 = generate_record(hdr,"chr1\t6700131\t.\tC\tCAAA\t0\tPASS\tMQ=60\tGT:GQ:GQX:DPI:AD:ADF:ADR:FT:PL\t0/0:93:93:41:37,1:18,0:19,1:PASS:0,96,545");
    std::deque<bcf1_t *> q;
    q.push_back(rec1);
    q.push_back(rec2);
    q.push_back(rec3);
    pair<std::deque<bcf1_t *>::iterator,std::deque<bcf1_t *>::iterator> i(q.begin(),q.end());
    auto collapsed_record = CollapseRecords(hdr,i);
 //   ggutils::print_variant(hdr,collapsed_record);   
    int32_t *ptr=nullptr,num_values=0;
    ASSERT_EQ(collapsed_record->n_allele,4);
    ASSERT_EQ(bcf_get_format_int32(hdr,collapsed_record,"AD",&ptr,&num_values),4);
    ASSERT_EQ(bcf_get_format_int32(hdr,collapsed_record,"PL",&ptr,&num_values),10);
    ASSERT_EQ(ptr[ 0 ], 189 );
    ASSERT_EQ(ptr[ 1 ], 0 );
    ASSERT_EQ(ptr[ 2 ], 197 );
    ASSERT_EQ(ptr[ 3 ], 297 );
    ASSERT_EQ(ptr[ 4 ], 108 );
    ASSERT_EQ(ptr[ 5 ], 750 );
    ASSERT_EQ(ptr[ 6 ], 285 );
    ASSERT_EQ(ptr[ 7 ], 96 );
    ASSERT_EQ(ptr[ 8 ], 393 );
    ASSERT_EQ(ptr[ 9 ], 734 );
}