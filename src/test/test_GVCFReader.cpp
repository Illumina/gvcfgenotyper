#include "gtest/gtest.h"
#include "GVCFReader.hpp"
#include "common.hpp"

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
    while(line!=NULL)
    {
	ofs << bcf_hdr_id2name(hdr,line->rid)<<":"<<line->pos+1<<":"<<line->d.allele[0]<<":"<<line->d.allele[1]<<std::endl;
	DepthBlock db;
	g.get_depth(line->rid,line->pos,line->pos+1,db);
	std::cerr << "depth at " << line->pos+1<< " = " << db._dp  << std::endl;	
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
