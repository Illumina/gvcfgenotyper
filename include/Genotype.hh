//
// Created by O'Connell, Jared on 11/2/17.
//

#ifndef GVCFGENOTYPER_GENOTYPE_HH
#define GVCFGENOTYPER_GENOTYPE_HH

extern "C" {
#include <htslib/vcf.h>
}

#include "utils.hh"
#include "multiAllele.hh"
#include  <deque>
//#include <pair>

class Genotype
{
public:
    Genotype(bcf_hdr_t  *header, bcf1_t *record);
    Genotype(int ploidy, int num_allele);
    Genotype marginalise(const int index);
    Genotype(bcf_hdr_t *sample_header,pair<std::deque<bcf1_t *>::iterator,std::deque<bcf1_t *>::iterator> & sample_variants,multiAllele & alleles_to_map);
    ~Genotype();
    void setDepthFromAD();
    void print();//prints information to stderr. mainly for debugging.
    int update_bcf1_t(int sample_index, bcf_hdr_t *header, bcf1_t *record);
    int propagate_format_fields(int sample_index,int ploidy,int *gt,int *gq,int *gqx,int *dp,int *dpf,int *ad,int *adf,int *adr,int *pl);
    int get_gq();
    int get_gqx();
    int get_dp();
    int get_dpf();
    int get_gt(int index);
    int get_ploidy();
    int get_ad(int index);
    int get_adr(int index);
    int get_adf(int index);
    int get_pl(int g0,int g1);
    void set_dp_missing();
    bool is_dp_missing();
    int *_gt, *_ad, *_gq, *_dp, *_dpf, *_pl, *_adf, *_adr, *_gqx;
    int _num_allele, _num_pl, _ploidy, _num_gt, _num_ad, _num_adf, _num_adr,  _num_gq, _num_gqx, _num_dp, _num_dpf, _num_gl;
    std::vector<float> _gl;
    bool _has_pl, _adf_found, _adr_found;
    void PLfromGL();
    int get_mq();
    float get_qual();
private:
    void allocate(int ploidy, int num_allele);
    float _qual;
    int32_t _mq;
};

#endif //GVCFGENOTYPER_GENOTYPE_HH
