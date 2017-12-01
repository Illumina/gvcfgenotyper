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

//Genotype stores FORMAT/INFO fields from a VCF record for a single sample.
//It contains a number of helper functions to manipulate these format/info
//values and propagate them back to new (possible multiple sample) VCF record.
class Genotype
{
public:
    //Constructs a Genotype with values taken from record.
    Genotype(bcf_hdr_t const *header, bcf1_t *record);

    //Constructs an empty Genotype with memory allocated according the ploidy/num_allele.
    Genotype(int ploidy, int num_allele);

    //Constructs a Genotype with alleles from alleles_to_map and format/info values taken from sample_variants (which is a subset of alleles_to_map).
    Genotype(bcf_hdr_t *sample_header,pair<std::deque<bcf1_t *>::iterator,std::deque<bcf1_t *>::iterator> & sample_variants,multiAllele & alleles_to_map);
    ~Genotype();

    //Marginalises over all alleles != index, creating a symbolic X allele which holds all the depth and probability mass for marginalised out alleles. This is being deprecated.
    Genotype marginalise(int index);

    //Removes alleles in indices, adding their AD and PL values to the REF values.
    Genotype collapse_alleles_into_ref(vector<int> & indices);

    //Copies FORMAT fields from this object into a vcf_data_t object, ploidy is the ploidy of the vcf_data_t format.
    int propagate_format_fields(size_t sample_index,size_t ploidy,ggutils::vcf_data_t *format);

    //Copies FORMAT/INFO fields from this object into record.
    //Copies FORMAT/INFO fields from this object into record.
    int update_bcf1_t(bcf_hdr_t *header, bcf1_t *record);

    //Updates FORMAT/DP by summing FORMAT/AD. This is be cause FORMAT/DP is not assigned at indels.
    void setDepthFromAD();

    //This phreds the _gl values and copies them to _pl.
    void PLfromGL();

    //Zeroes DP and AD* values.
    void set_depth_to_zero();

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
    bool is_dp_missing();
    int *_gt, *_ad, *_gq, *_dp, *_dpf, *_pl, *_adf, *_adr, *_gqx;
    int _num_allele, _num_pl, _ploidy, _num_gt, _num_ad, _num_adf, _num_adr,  _num_gq, _num_gqx, _num_dp, _num_dpf, _num_gl;
    std::vector<float> _gl;
private:
    //Assigns memory according to ploidy/num_allele.
    void allocate(int ploidy, int num_allele);
    float _qual;
    int32_t _mq;
    bool _has_pl, _adf_found, _adr_found;
};

#endif //GVCFGENOTYPER_GENOTYPE_HH
