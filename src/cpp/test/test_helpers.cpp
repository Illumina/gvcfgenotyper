//
// Created by O'Connell, Jared on 11/1/17.
//
#include "test_helpers.hh"

bcf1_t *generate_record(bcf_hdr_t *hdr, kstring_t & vcfrow)
{
    bcf1_t *ret = bcf_init1();
    vcf_parse(&vcfrow,hdr,ret);
    bcf_unpack(ret,BCF_UN_ALL);
    free(vcfrow.s);
    return(ret);
}

bcf1_t *generate_record(bcf_hdr_t *hdr,const std::string & vcfrow)
{
    kstring_t str = {0,0,nullptr};
    kputs(vcfrow.c_str(),&str);
    return generate_record(hdr,str);
}

// generates a record with a GT
bcf1_t *generate_record(bcf_hdr_t *hdr,int rid,int pos,const std::string & alleles, const std::string& gt)
{
    kstring_t str = {0,0,nullptr};
    size_t comma = alleles.find(",");
    assert(comma<alleles.size());
    std::string ref = alleles.substr(0,comma);
    std::string alt = alleles.substr(comma+1,alleles.size());
    ksprintf(&str,"%s\t%d\t.\t%s\t%s\t1691\tPASS\t.\tGT\t%s",bcf_hdr_int2id(hdr,BCF_DT_CTG,rid),pos,ref.c_str(),alt.c_str(),gt.c_str());
    return generate_record(hdr,str);
}

bcf1_t *generate_record(bcf_hdr_t *hdr,int rid,int pos,const std::string & alleles)
{
    return generate_record(hdr,rid,pos,alleles,"0/1");
}

bcf_hdr_t *get_header()
{
    //FIXME: we should store the header as a big string within the source code
    std::string gvcf_file_name = g_testenv->getBasePath() + "/../test/test2/NA12887_S1.vcf.gz";
    auto bcf_header=bcf_hdr_read(hts_open(gvcf_file_name.c_str(), "r"));
    bcf_hdr_append(bcf_header, "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Sample filter, 'PASS' indicates that all single sample filters passed for this sample\">");
    bcf_hdr_sync(bcf_header);
    return(bcf_header);
}

void update_record(bcf_hdr_t *hdr,int rid,int pos,const std::string & alleles,bcf1_t *record)
{
    bcf_clear1(record);
    kstring_t str = {0,0,nullptr};
    size_t comma = alleles.find(",");
    assert(comma<alleles.size());
    std::string ref = alleles.substr(0,comma);
    std::string alt = alleles.substr(comma+1,alleles.size());
    ksprintf(&str,"%s\t%d\t.\t%s\t%s\t1691\tPASS\t.\tGT\t0/1",bcf_hdr_int2id(hdr,BCF_DT_CTG,rid),pos,ref.c_str(),alt.c_str());
    vcf_parse(&str,hdr,record);
    bcf_unpack(record,BCF_UN_ALL);
    free(str.s);
}


