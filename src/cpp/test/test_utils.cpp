#include <htslib/vcf.h>
#include "test_helpers.hh"

#include "ggutils.hh"

TEST(UtilTest, comparators)
{
    std::string gvcf_file_name = g_testenv->getBasePath() + "/../test/NA12877.tiny.vcf.gz";

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

    update_record(hdr,1,100,"TA,TAA",record1);
    update_record(hdr,1,100,"T,TA" ,record2);
    ASSERT_TRUE( ggutils::bcf1_equal(record1, record2));

    update_record(hdr,1,100,"TAA,TA",record1);
    update_record(hdr,1,100,"TA,T" ,record2);
    ASSERT_TRUE( ggutils::bcf1_equal(record1, record2));

    update_record(hdr,1,100,"TAA,TA,T",record1);
    update_record(hdr,1,100,"TA,T" ,record2);
    ASSERT_TRUE( ggutils::bcf1_equal(record1, record2));

    update_record(hdr,1,100,"TAA,T,TA",record1);
    update_record(hdr,1,100,"TA,T" ,record2);
    ASSERT_FALSE( ggutils::bcf1_equal(record1, record2));

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
    ASSERT_EQ(105,ggutils::get_gl_index(0,14));
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

TEST(UtilTest,bcf1AlleleSwapDiploid)
{
    auto hdr = get_header();
    auto record1 = generate_record(hdr, "chr1\t5420\t.\tC\tA,T,G\t100\tPASS\tMQ=50;AF1000G=.13\tGT:GQ:DP:DPF:AD:PL:SB\t1/3:30:16:0:0,12,0,4:0,1,2,3,4,5,6,7,8,9:1.01");
    int32_t *pl=nullptr,*pl_new=nullptr,num_pl=0,num_pl_new=0;
    bcf_get_format_int32(hdr,record1,"PL",&pl,&num_pl);
    ggutils::print_variant(hdr,record1);
    for(int i=1;i<record1->n_allele;i++)
    {
        bcf1_t *record2=bcf_dup(record1);
        bcf_unpack(record2,BCF_UN_ALL);
        ggutils::bcf1_allele_swap(hdr,record2,i,1);
        ggutils::print_variant(hdr,record2);
        bcf_get_format_int32(hdr,record2,"PL",&pl_new,&num_pl_new);
        ASSERT_EQ(pl[ggutils::get_gl_index(i,0)],pl_new[ggutils::get_gl_index(1,0)]);
        ASSERT_EQ(pl[ggutils::get_gl_index(i,1)],pl_new[ggutils::get_gl_index(1,i)]);
        ASSERT_EQ(pl[ggutils::get_gl_index(i,i)],pl_new[ggutils::get_gl_index(1,1)]);
        for(int j=2;j<record1->n_allele;j++)
        {
            if(i!=j) {ASSERT_EQ(pl[ggutils::get_gl_index(i,j)],pl_new[ggutils::get_gl_index(1,j)]);}
        }
        bcf_destroy1(record2);
    }
    free(pl);
    free(pl_new);
}

TEST(UtilTest,bcf1AlleleSwapHaploid)
{
    auto hdr = get_header();
    auto record1 = generate_record(hdr,"chr1\t1\t.\tA\tC,G\t60\tSiteConflict\tMQ=51\tGT:GQ:GQX:DP:DPF:AD:ADF:ADR:SB:FT:PL\t2:17:8:6:10:0,3,3:0,0,2:0,3,1:-7.3:SiteConflict:95,17,0");

    int32_t *pl=nullptr,*pl_new=nullptr,num_pl=0,num_pl_new=0;
    bcf_get_format_int32(hdr,record1,"PL",&pl,&num_pl);
    ggutils::print_variant(hdr,record1);
    bcf1_t *record2=bcf_dup(record1);
    bcf_unpack(record2,BCF_UN_ALL);
    ggutils::bcf1_allele_swap(hdr,record2,2,1);
    ggutils::print_variant(hdr,record2);
    bcf_get_format_int32(hdr,record2,"PL",&pl_new,&num_pl_new);
    ASSERT_EQ(pl[0],pl_new[0]);
    ASSERT_EQ(pl[1],pl_new[2]);
    ASSERT_EQ(pl[2],pl_new[1]);
    free(pl);
    free(pl_new);
}

TEST(UtilTest,rightTrim)
{
    std::string ref="TA";
    std::string alt="TAA";
    size_t a,b;
    ggutils::right_trim(ref.c_str(),alt.c_str(),a,b);
    ASSERT_EQ(a,(size_t)1);
    ASSERT_EQ(b,(size_t)2);

    ref="T";
    alt="TA";
    ggutils::right_trim(ref.c_str(),alt.c_str(),a,b);
    ASSERT_EQ(a,(size_t)1);
    ASSERT_EQ(b,(size_t)2);

    ref="TTTAAATTCT";
    alt="TTTAAATACT";
    ggutils::right_trim(ref.c_str(),alt.c_str(),a,b);
    ASSERT_EQ(a,ref.size()-2);
    ASSERT_EQ(b,alt.size()-2);

    ref="T";
    alt="A";
    ggutils::right_trim(ref.c_str(),alt.c_str(),a,b);
    ASSERT_EQ(a,(size_t)1);
    ASSERT_EQ(b,(size_t)1);

    ref="T";
    alt="T";
    ggutils::right_trim(ref.c_str(),alt.c_str(),a,b);
    ASSERT_EQ(a,(size_t)1);
    ASSERT_EQ(b,(size_t)1);
}

TEST(UtilTest,findAllele)
{
    auto hdr = get_header();
    auto target = generate_record(hdr, "chr1\t1\t.\tC\tA,T,G\t0\t.\t.\tGT\t1/3");
    auto record2 = generate_record(hdr, "chr1\t1\t.\tC\tA,T,G\t0\t.\t.\tGT\t1/3");
    for(int i=1;i<target->n_allele;i++)
        ASSERT_EQ(ggutils::find_allele(target,record2,i),i);

    record2 = generate_record(hdr, "chr1\t1\t.\tC\tG\t0\t.\t.\tGT\t0/1");
    ASSERT_EQ(ggutils::find_allele(target,record2,1),3);

    auto target2 = generate_record(hdr,"chr1\t1\t.\tAT\tATTT,ATT,A\t589\tPASS\tMQ=60\tGT\t1/1");
    auto query1 = generate_record(hdr,"chr1\t1\t.\tAT\tATT,ATTT,A\t589\tPASS\tMQ=60\tGT\t1/1");
    auto query2 = generate_record(hdr,"chr1\t1\t.\tAT\tA,ATTT,ATT\t589\tPASS\tMQ=60\tGT\t1/1");

    for(int i=1;i<target->n_allele;i++)
        ASSERT_EQ(ggutils::find_allele(target2,target2,i),i);

    ASSERT_EQ(ggutils::find_allele(target2,query1,1),2);
    ASSERT_EQ(ggutils::find_allele(target2,query1,2),1);
    ASSERT_EQ(ggutils::find_allele(target2,query1,3),3);

    ASSERT_EQ(ggutils::find_allele(target2,query2,1),3);
    ASSERT_EQ(ggutils::find_allele(target2,query2,2),1);
    ASSERT_EQ(ggutils::find_allele(target2,query2,3),2);
}

TEST(UtilTest,collapseLikelihoods1)
{
    int p[][3] = { {189,0,197},{0,108,561},{0,96,545} };
    int num_allele = 4;
    int num_gl = ggutils::get_number_of_likelihoods(2,num_allele);
    std::vector< std::vector<int> > gl;
    for(int i=1;i<4;i++)
    {
        std::vector<int> tmp(num_gl,bcf_int32_missing);
        tmp[ggutils::get_gl_index(0,0)] = p[i-1][0];
        tmp[ggutils::get_gl_index(0,i)] = p[i-1][1];
        tmp[ggutils::get_gl_index(i,i)] = p[i-1][2];
        gl.push_back(tmp);
    }
    std::vector<int> new_gl;
    ggutils::collapse_gls(2,num_allele,gl,new_gl);
    ASSERT_EQ(new_gl[ 0 ], 189 );
    ASSERT_EQ(new_gl[ 1 ], 0 );
    ASSERT_EQ(new_gl[ 2 ], 197 );
    ASSERT_EQ(new_gl[ 3 ], 297 );
    ASSERT_EQ(new_gl[ 4 ], 108 );
    ASSERT_EQ(new_gl[ 5 ], 750 );
    ASSERT_EQ(new_gl[ 6 ], 285 );
    ASSERT_EQ(new_gl[ 7 ], 96 );
    ASSERT_EQ(new_gl[ 8 ], 393 );
    ASSERT_EQ(new_gl[ 9 ], 734 );
}

TEST(UtilTest,collapseLikelihoods2)
{
    int p[][10] = { {189,0,197,bcf_int32_missing,bcf_int32_missing,bcf_int32_missing,bcf_int32_missing,bcf_int32_missing,bcf_int32_missing,bcf_int32_missing},
                   {0,bcf_int32_missing,bcf_int32_missing,108,bcf_int32_missing,561,96,bcf_int32_missing,bcf_int32_missing,545 } };

    int num_allele = 4;
    int num_gl = ggutils::get_number_of_likelihoods(2,num_allele);
    std::vector< std::vector<int> > gl;
    gl.emplace_back(p[0],p[0]+num_gl);
    gl.emplace_back(p[1],p[1]+num_gl);
    std::vector<int> new_gl;
    ggutils::collapse_gls(2,num_allele,gl,new_gl);
//    for(auto it=new_gl.begin();it!=new_gl.end();it++) std::cerr << *it <<"\t";    std::cerr<<std::endl;
    ASSERT_EQ(new_gl[ 0 ], 189 );
    ASSERT_EQ(new_gl[ 1 ], 0 );
    ASSERT_EQ(new_gl[ 2 ], 197 );
    ASSERT_EQ(new_gl[ 3 ], 297 );
    ASSERT_EQ(new_gl[ 4 ], 108 );
    ASSERT_EQ(new_gl[ 5 ], 750 );
    ASSERT_EQ(new_gl[ 6 ], 285 );
    ASSERT_EQ(new_gl[ 7 ], 96 );
    ASSERT_EQ(new_gl[ 8 ], 297 );
    ASSERT_EQ(new_gl[ 9 ], 734 );
}

TEST(UtilTest,addAllele1)
{
    auto hdr = get_header();
    auto v1 = generate_record(hdr, "chr1\t1\t.\tC\tA\t0\t.\t.\tGT\t0/1");
    auto v2 = generate_record(hdr, "chr1\t1\t.\tC\tT\t0\t.\t.\tGT\t0/1");
    ggutils::add_allele(hdr,v1,v2,1);
    ASSERT_STREQ(v1->d.allele[0],"C");
    ASSERT_STREQ(v1->d.allele[1],"A");
    ASSERT_STREQ(v1->d.allele[2],"T");
}

TEST(UtilTest,addAllele2)
{
    auto hdr = get_header();
    auto v1 = generate_record(hdr, "chr1\t1\t.\tCA\tC\t0\t.\t.\tGT\t0/1");
    auto v2 = generate_record(hdr, "chr1\t1\t.\tCAA\tC\t0\t.\t.\tGT\t0/1");
    ggutils::add_allele(hdr,v1,v2,1);
    ASSERT_STREQ(v1->d.allele[0],"CAA");
    ASSERT_STREQ(v1->d.allele[1],"CA");
    ASSERT_STREQ(v1->d.allele[2],"C");
}

TEST(UtilTest,addAllele3)
{
    auto hdr = get_header();
    auto v1 = generate_record(hdr, "chr1\t1\t.\tCAA\tC\t0\t.\t.\tGT\t0/1");
    auto v2 = generate_record(hdr, "chr1\t1\t.\tCA\tC\t0\t.\t.\tGT\t0/1");
    ggutils::add_allele(hdr,v1,v2,1);
    ASSERT_STREQ(v1->d.allele[0],"CAA");
    ASSERT_STREQ(v1->d.allele[1],"C");
    ASSERT_STREQ(v1->d.allele[2],"CA");
}

TEST(UtilTest,addAllele4)
{
    auto hdr = get_header();
    auto v1 = generate_record(hdr,"chr21\t9437597\t.\tACC\tCTCCCCGCCGCCGTGGCTTTTTGACA,CTCCCCGCCGCCGTGGCTTTTTGACACCGCCGCCGCGGCTTTTGGTCC\t42\tPASS\tCIGAR=1M3D26I,1M1D46I2M;RU=.,.;REFREP=.,.;IDREP=.,.\tGT:GQ:GQX:DPI:AD\t1/2:93:53:21:12,3,3");
    auto v2 = generate_record(hdr,"chr21\t9437597\t.\tA\tC\t0\tSiteConflict;LowGQX;HighDPFRatio\t.\tGT:GQX:DP:DPF:AD\t0/1:25:4:8:3,1");
    ggutils::add_allele(hdr,v2,v1,2);
    ASSERT_STREQ(v2->d.allele[2],"CTCCCCGCCGCCGTGGCTTTTTGACACCGCCGCCGCGGCTTTTGGT");
    ggutils::add_allele(hdr,v2,v1,1);
    ggutils::add_allele(hdr,v2,v2,1);
    ggutils::add_allele(hdr,v2,v1,1);
    ggutils::add_allele(hdr,v2,v1,2);

///    ggutils::print_variant(hdr,v2);
    ASSERT_STREQ(v2->d.allele[0],"ACC");
    ASSERT_STREQ(v2->d.allele[1],"CCC");
    ASSERT_STREQ(v2->d.allele[2],"CTCCCCGCCGCCGTGGCTTTTTGACACCGCCGCCGCGGCTTTTGGTCC");
    ASSERT_STREQ(v2->d.allele[3],"CTCCCCGCCGCCGTGGCTTTTTGACA");

}

TEST(UtilTest,median)
{
    int x[5] = {10,5,3,2,289};
    ASSERT_FLOAT_EQ(5.,ggutils::median(x,5));
    int y[6] = {10,5,3,2,289,1000000};
    ASSERT_FLOAT_EQ(7.5,ggutils::median(y,6));
}

TEST(UtilTest,fisherSB)
{
    std::vector<float> p;
    int adf[3] = {10,25,1};
    int adr[3] = {13,26,40};
    ggutils::fisher_sb_test(adf,adr,3,p);
    ASSERT_FLOAT_EQ(p[1],4.187761);
    ASSERT_FLOAT_EQ(p[0],0.09571717);
}