//
// Created by O'Connell, Jared on 10/23/17.
//

#ifndef GVCFGENOTYPER_TEST_UTILD_H
#define GVCFGENOTYPER_TEST_UTILD_H



static bcf1_t *generate_record(bcf_hdr_t *hdr,int rid,int pos,const string & alleles)
{
    bcf1_t *ret = bcf_init1();
    kstring_t str = {0,0,nullptr};
    int comma = alleles.find(",");
    assert(comma<alleles.size());
    string ref = alleles.substr(0,comma);
    string alt = alleles.substr(comma+1,alleles.size());
    ksprintf(&str,"%s\t%d\t.\t%s\t%s\t1691\tPASS\t.\tGT\t0/1",bcf_hdr_int2id(hdr,BCF_DT_CTG,rid),pos,ref.c_str(),alt.c_str());
    vcf_parse(&str,hdr,ret);
    bcf_unpack(ret,BCF_UN_ALL);
    free(str.s);
    return(ret);
}


static bcf1_t *generate_record(bcf_hdr_t *hdr,const string & vcfrow)
{
    bcf1_t *ret = bcf_init1();
    kstring_t str = {0,0,nullptr};
    kputs(vcfrow.c_str(),&str);
    vcf_parse(&str,hdr,ret);
    bcf_unpack(ret,BCF_UN_ALL);
    free(str.s);
    return(ret);
}


static bcf_hdr_t *get_header()
{
    std::string gvcf_file_name = g_testenv->getBasePath() + "/data/test2/NA12887_S1.vcf.gz";
    return(bcf_hdr_read(hts_open(gvcf_file_name.c_str(), "r")));
}

static void update_record(bcf_hdr_t *hdr,int rid,int pos,const string & alleles,bcf1_t *record)
{
    bcf_clear1(record);
    kstring_t str = {0,0,nullptr};
    int comma = alleles.find(",");
    assert(comma<alleles.size());
    string ref = alleles.substr(0,comma);
    string alt = alleles.substr(comma+1,alleles.size());
    ksprintf(&str,"%s\t%d\t.\t%s\t%s\t1691\tPASS\t.\tGT\t0/1",bcf_hdr_int2id(hdr,BCF_DT_CTG,rid),pos,ref.c_str(),alt.c_str());
    vcf_parse(&str,hdr,record);
    bcf_unpack(record,BCF_UN_ALL);
    free(str.s);
}


#endif //GVCFGENOTYPER_TEST_UTILD_H
