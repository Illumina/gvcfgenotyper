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
#include "spdlog.h"

//Genotype stores FORMAT/INFO fields from a VCF record for a single sample.
//It contains a number of helper functions to manipulate these format/info
//values and propagate them back to new (possible multiple sample) VCF record.
class Genotype
{
public:
    //Constructs a Genotype with values taken from record.
    Genotype(bcf_hdr_t const *header, bcf1_t *record);
    void Init(bcf_hdr_t const *header, bcf1_t *record);

    //Constructs an empty Genotype with memory allocated according the ploidy/num_allele.
    Genotype(int ploidy, int num_allele);

    //Constructs a Genotype with alleles from alleles_to_map and format/info values taken from sample_variants (which contains subset of alleles_to_map).
    //Handle "conflicts" where allele and genotype combinations conflict with one another in a rudimentary but sane way.
    Genotype(bcf_hdr_t *sample_header,bcf1_t *sample_variants,multiAllele & alleles_to_map);

    Genotype(bcf_hdr_t *sample_header,pair<std::deque<bcf1_t *>::iterator,std::deque<bcf1_t *>::iterator> & sample_variants);

    ~Genotype();

    //Removes alleles in indices, adding their AD and PL values to the REF values.
    void CollapseAllelesIntoRef(vector<int> &indices, Genotype &out);

    //Copies FORMAT fields from this object into a vcf_data_t object, ploidy is the ploidy of the vcf_data_t format.
    int PropagateFormatFields(size_t sample_index, size_t ploidy, ggutils::vcf_data_t *format);

    //Copies FORMAT/INFO fields from this object into record.
    int UpdateBcfRecord(bcf_hdr_t *header, bcf1_t *record);

    //Updates FORMAT/DP by summing FORMAT/AD. This is be cause FORMAT/DP is not assigned at indels.
    void SetDepthFromAd();

    //These phreds the _gl values and copies them to _pl and vice versa.
    void SetPlFromGl();
    void SetGlFromPl();

    //Zeroes DP and AD* values.
    void SetDepthToZero();
    void SetGtToHomRef();

    //takes a haploid call and makes it diploid. For genotypes, 0 -> 0/0 and 1 -> 1/1 etc.
    void MakeDiploid();
    //sets FORMAT/GT to argmax(GL)
    void CallGenotype();
    bool IsDpMissing();
    bool IsGtMissing() {return bcf_gt_is_missing(_gt[0]);}

    //accessors/mutators
    int gq();
    int gqx();
    int dp();
    int dpf();
    int ad(int index);
    int adr(int index);
    int adf(int index);
    int gt(int index);
    const char *filter();
    
    int mq();
    int ploidy();
    int num_allele();
    int num_pl() {return _num_pl;};
    int pl(int g0, int g1);
    int pl(int g0);
    float gl(int g0, int g1);
    float gl(int g0);
    float qual();

    void SetDp(int val);
    void SetDpf(int val);
    void SetGq(int val);
    void SetGqx(int val);
    void SetAd(int val,int index);
    void SetAdf(int val,int index);
    void SetAdr(int val,int index);
    void SetPl(std::vector<int> & val);
    bool HasPl() {return _has_pl;};
    bool HasAdf() {return _adf_found;};
    bool HasAdr() {return _adr_found;};

    int32_t *_gt=nullptr, *_ad=nullptr, *_gq=nullptr, *_dp=nullptr, *_dpf=nullptr, *_pl=nullptr, *_adf=nullptr, *_adr=nullptr, *_gqx=nullptr;
    int _num_allele, _ploidy, _num_pl=0, _num_gt=0, _num_ad=0, _num_adf=0, _num_adr=0,  _num_gq=0, _num_gqx=0, _num_dp=0, _num_dpf=0,_num_filter=0;
    std::vector<float> _gl;
    char *_filter=nullptr;//stores the FILTER column as a string.
    
private:
    //Assigns memory according to ploidy/num_allele.
    void allocate(int ploidy, int num_allele);
    void init_logger();
    float _qual;
    int32_t _mq;
    bool _has_pl, _adf_found, _adr_found;
    std::shared_ptr<spdlog::logger> _lg;
};

#endif //GVCFGENOTYPER_GENOTYPE_HH
