//
// Created by O'Connell, Jared on 11/2/17.
//

#ifndef GVCFGENOTYPER_GENOTYPE_HH
#define GVCFGENOTYPER_GENOTYPE_HH

#include <utility>
#include <deque>


extern "C" {
#include <htslib/vcf.h>
}

#include "multiAllele.hh"
#include "ggutils.hh"


class Genotype
{
public:
    Genotype(bcf_hdr_t const *header, bcf1_t *record);
    Genotype(int ploidy, int num_allele);
    Genotype(bcf_hdr_t *sample_header,pair<std::deque<bcf1_t *>::iterator,std::deque<bcf1_t *>::iterator> & sample_variants,multiAllele & alleles_to_map);
    int propagate_format_fields(size_t sample_index,size_t ploidy,ggutils::vcf_data_t *format);
    Genotype marginalise(int index);
    Genotype collapse_alleles_into_ref(vector<int> & indices);
    ~Genotype();
    void setDepthFromAD();
    int update_bcf1_t(bcf_hdr_t *header, bcf1_t *record);
    int get_gq();
    int get_gqx();
    int get_dp();
    int get_dpf();
    int get_ad(int index);
    int get_adr(int index);
    int get_adf(int index);
    int get_gt(int index);
    int get_mq();
    int get_ploidy();
    int get_pl(int g0,int g1);
    int get_pl(int g0);
    float get_qual();
    void PLfromGL();
    void set_dp_missing();
    bool is_dp_missing();
    int *_gt, *_ad, *_gq, *_dp, *_dpf, *_pl, *_adf, *_adr, *_gqx;
    int _num_allele, _num_pl, _ploidy, _num_gt, _num_ad, _num_adf, _num_adr,  _num_gq, _num_gqx, _num_dp, _num_dpf, _num_gl;
    std::vector<float> _gl;
    bool _has_pl, _adf_found, _adr_found;
private:
    void allocate(int ploidy, int num_allele);
    float _qual;
    int32_t _mq;
};

#endif //GVCFGENOTYPER_GENOTYPE_HH
