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
    Genotype(bcf_hdr_t *sample_header,pair<std::deque<bcf1_t *>::iterator,std::deque<bcf1_t *>::iterator> & sample_variants,
             multiAllele & alleles_to_map);
    ~Genotype();

    //Removes alleles in indices, adding their AD and PL values to the REF values.
    void collapse_alleles_into_ref(vector<int> & indices,Genotype & out);

    //Copies FORMAT fields from this object into a vcf_data_t object, ploidy is the ploidy of the vcf_data_t format.
    int propagate_format_fields(size_t sample_index,size_t ploidy,ggutils::vcf_data_t *format);

    //Copies FORMAT/INFO fields from this object into record.
    //Copies FORMAT/INFO fields from this object into record.
    int update_bcf1_t(bcf_hdr_t *header, bcf1_t *record);

    //Updates FORMAT/DP by summing FORMAT/AD. This is be cause FORMAT/DP is not assigned at indels.
    void set_depth_from_ad();

    //These phreds the _gl values and copies them to _pl and vice versa.
    void PLfromGL();
    void GLfromPL();

    //Zeroes DP and AD* values.
    void set_depth_to_zero();
    void set_gt_to_homref();

    //takes a haploid call and makes it diploid. For genotypes, 0 -> 0/0 and 1 -> 1/1 etc.
    void make_diploid();

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
    int get_num_allele();
    int get_pl(int g0,int g1);
    int get_pl(int g0);
    float get_qual();

    void set_dp(int val);
    void set_dpf(int val);
    void set_gq(int val);
    void set_gqx(int val);

    bool is_dp_missing();
    int32_t *_gt=nullptr, *_ad=nullptr, *_gq=nullptr, *_dp=nullptr, *_dpf=nullptr, *_pl=nullptr, *_adf=nullptr, *_adr=nullptr, *_gqx=nullptr;
    int _num_allele, _ploidy, _num_pl=0, _num_gt=0, _num_ad=0, _num_adf=0, _num_adr=0,  _num_gq=0, _num_gqx=0, _num_dp=0, _num_dpf=0;
    std::vector<float> _gl;
private:
    //Assigns memory according to ploidy/num_allele.
    void allocate(int ploidy, int num_allele);
    float _qual;
    int32_t _mq;
    bool _has_pl, _adf_found, _adr_found;
};

#endif //GVCFGENOTYPER_GENOTYPE_HH
