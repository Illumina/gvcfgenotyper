//
// Created by O'Connell, Jared on 10/23/17.
//

#ifndef GVCFGENOTYPER_TEST_UTILD_H
#define GVCFGENOTYPER_TEST_UTILD_H


static bcf1_t *generate_record(bcf_hdr_t *hdr,int rid,int pos,const string & alleles)
{
    bcf1_t *ret = bcf_init1();
    ret->rid = rid;
    ret->pos = pos;
    bcf_update_alleles_str(hdr, ret, alleles.c_str());
    return(ret);
}

static bcf_hdr_t *get_header()
{
    std::string gvcf_file_name = g_testenv->getBasePath() + "/data/NA12877.tiny.vcf.gz";
    return(bcf_hdr_read(hts_open(gvcf_file_name.c_str(), "r")));
}
#endif //GVCFGENOTYPER_TEST_UTILD_H
