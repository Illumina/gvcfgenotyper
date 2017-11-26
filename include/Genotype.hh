//
// Created by O'Connell, Jared on 11/2/17.
//

#ifndef GVCFGENOTYPER_GENOTYPE_HH
#define GVCFGENOTYPER_GENOTYPE_HH

extern "C" {
#include <htslib/vcf.h>
}

#include "utils.hh"


class Genotype
{
public:
    Genotype(bcf_hdr_t const *header, bcf1_t *record);
    Genotype(int ploidy, int num_allele);
    Genotype marginalise(int index);
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
    void set_dp_missing();
    bool is_dp_missing();
    int *_gt, *_ad, *_gq, *_dp, *_dpf, *_pl, *_adf, *_adr, *_gqx;
    int _num_allele, _num_pl, _ploidy, _num_gt, _num_ad, _num_adf, _num_adr,  _num_gq, _num_gqx, _num_dp, _num_dpf, _num_gl;
    std::vector<float> _gl;
    bool _has_pl, _adf_found, _adr_found;
};

#endif //GVCFGENOTYPER_GENOTYPE_HH
