#include "gtest/gtest.h"
#include "GVCFReader.hpp"
#include "common.hpp"


TEST(DepthBlock,depth_block_unit_tests)
{
    DepthBlock db1(0,100,199,20,1,30);
    DepthBlock db2(0,100,100,20,1,30);
    ASSERT_EQ(db1.intersect_size(db2),1);
    ASSERT_EQ(db1.intersect_size(db2),db2.intersect_size(db1));

    DepthBlock db3(0,300,401,20,1,30);
    ASSERT_EQ(db1.intersect_size(db3),0);

    DepthBlock db4(1,300,401,20,1,30);
    ASSERT_EQ(db1.intersect_size(db4),0);

    DepthBlock db5(0,190,401,20,1,30);
    ASSERT_EQ(db1.intersect_size(db5),10);

    DepthBlock db6(0,100,100,20,1,30);
    DepthBlock db7(0,100,100,20,1,30);
    ASSERT_EQ(db6.intersect_size(db7),1);
}

//GVCFReader should match this VID output
//bcftools norm -m -any data/NA12877.tiny.vcf.gz | bcftools norm -f data/tiny.ref.fa | bcftools query -i 'ALT!="."' -f '%CHROM:%POS:%REF:%ALT\n' > data/NA12877.tiny.vcf.gz.expected 
TEST(GVCFReader,tiny_gvcf_example)
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
    int buffer_size=200;
    GVCFReader g(gvcf_file_name,ref_file_name,buffer_size);    
    const bcf_hdr_t *hdr = g.getHeader();    
    bcf1_t *line = g.pop();
    int32_t *dp=NULL,nval=0;
    DepthBlock db;
    while(line!=NULL)
    {
	ofs << bcf_hdr_id2name(hdr,line->rid)<<":"<<line->pos+1<<":"<<line->d.allele[0]<<":"<<line->d.allele[1]<<std::endl;
	if(is_snp(line))
	{
	    g.get_depth(line->rid,line->pos,line->pos,db);
	    if(bcf_get_format_int32(hdr,line, "DP", &dp,&nval)==1)
	    {
		ASSERT_EQ(db._dp,*dp);		
	    }
	}
	assert(line->n_allele==2);
	bcf_destroy(line);
	line = g.pop();
    }
    ofs.close();
    ASSERT_TRUE(g.empty());
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


