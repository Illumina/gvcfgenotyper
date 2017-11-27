//
// Created by O'Connell, Jared on 11/2/17.
//

#ifndef GVCFGENOTYPER_NORMALISER_HH
#define GVCFGENOTYPER_NORMALISER_HH

extern "C" {
#include <htslib/vcf.h>
#include "vcfnorm.h"
}

#include "ggutils.hh"
#include "Genotype.hh"

//vcfnorm stuff
#define ERR_DUP_ALLELE      -2
#define ERR_REF_MISMATCH    -1
#define ERR_OK              0
#define ERR_SYMBOLIC        1
#define CHECK_REF_WARN 1
#define MROWS_SPLIT 1
#define MROWS_MERGE  2
//end vcf norm stuff

int mnp_split(bcf1_t *record_to_split, bcf_hdr_t *header, vector<bcf1_t *> &output);

// This class is problematic, with several members that are pointers to other
// classes it should define copy ctor and assignment operator as well. Ideally std::unique_ptr as well.
//this basically wraps bcftools norm in a class.
//TODO: investigate replacing this with invariant components
class Normaliser
{
public:
    Normaliser(const std::string &ref_fname, bcf_hdr_t *hdr);

    ~Normaliser();

    void unarise(bcf1_t *rec, std::vector<bcf1_t *>& atomised_variants);
    void multi_split(bcf1_t *bcf_record_to_split,vector<bcf1_t*>& split_variants);

private:
    char _symbolic_allele[2];
    args_t *_norm_args;
    bcf_hdr_t *_hdr;
};


#endif //GVCFGENOTYPER_NORMALISER_HH
