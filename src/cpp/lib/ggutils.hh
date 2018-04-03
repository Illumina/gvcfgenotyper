#pragma once

#define __STDC_LIMIT_MACROS

#include <stdint.h>
#include "math.h"
#include <stdlib.h>
#include <ctime>

#include <limits>
#include <set>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <deque>


using namespace std;

extern "C" {
#include <htslib/bgzf.h>
#include <htslib/tbx.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <htslib/synced_bcf_reader.h>
#include "htslib/kfunc.h"
}


namespace ggutils
{
    //Simple struct to hold our default FORMAT fields.
    struct vcf_data_t
    {
        int32_t *pl,*ad,*adf,*adr,*gt,*gq,*gqx,*dp,*dpf,*ps;
        size_t ploidy,num_allele,num_sample,num_ad,num_pl;
        vcf_data_t(size_t ploidy,size_t num_allele,size_t num_sample);
        void resize(size_t num_alleles);
        void set_missing();
        ~vcf_data_t();
    };

    void init_vcf_data(size_t ploidy,size_t num_allele,size_t num_sample,vcf_data_t & record);
    void destroy_vcf_data(vcf_data_t & record);

    int read_text_file(const string &fname, vector<string> &output);

    string record2string(bcf_hdr_t const *header, bcf1_t *record);
    void print_variant(bcf_hdr_t const *header, bcf1_t *record);

    void print_variant(bcf1_t *record);

    int32_t *zeros(int n);

    int32_t *assign_bcf_int32_missing(int32_t *ptr,size_t n);

    bool fileexists(const string &fname);

    inline void die(const string &s)
    {
        cerr << "ERROR: " << s << "\nExiting..." << endl;
        exit(1);
    }

    int strsplit(const string &input, const char split, vector<string> &out);

    string join(const vector<string> &input, const string &delim);

    vector<int> match(const vector<string> &x, const vector<string> &y);

    inline void warn(const string &s)
    {
        cerr << "WARNING: " << s << endl;
    }

    int copy_contigs(const bcf_hdr_t *src, bcf_hdr_t *dst);

    //simple dumps text from fname into output
    int read_text_file(const string &fname, vector<string> &output);

    bool is_variant(bcf1_t const *record);

    bool is_snp(bcf1_t *record);

    bool is_deletion(bcf1_t *record);

    bool is_insertion(bcf1_t *record);

    bool is_snp(bcf1_t *record);


    bool is_complex(bcf1_t *record);

    int get_variant_rank(bcf1_t *record);

    int get_end_of_gvcf_block(bcf_hdr_t *header, bcf1_t *record);


    int get_ploidy(bcf_hdr_t *header, bcf1_t *record);

    bool is_hom_ref(const bcf_hdr_t* header, bcf1_t* record);

    int get_variant_rank(bcf1_t *record);


    int get_end_of_variant(bcf1_t *record);

    bool bcf1_equal(bcf1_t *a, bcf1_t *b); //checks if the first ALT allele is equivalent in a/b
    bool bcf1_all_equal(bcf1_t *a, bcf1_t *b); //checks if all ALT alleles are equivalent

    bool bcf1_less_than(bcf1_t *a, bcf1_t *b);

    bool bcf1_greater_than(bcf1_t *a, bcf1_t *b);

    bool bcf1_leq(bcf1_t *a, bcf1_t *b);

    bool bcf1_geq(bcf1_t *a, bcf1_t *b);

    bool bcf1_not_equal(bcf1_t *a, bcf1_t *b);

    size_t get_number_of_likelihoods(int ploidy, int num_allele);

    int phred(float l);

    float unphred(int pl);

    int choose(int n, int k);

    //these are simple htslib wrappers shortcuts to get scalar ints/floats (not appropriate for arrays)
    int bcf1_get_one_info_float(const bcf_hdr_t *header, bcf1_t *record, const char *tag,float &output);
    int bcf1_get_one_format_float(const bcf_hdr_t *header, bcf1_t *record, const char *tag,float &output);
    int bcf1_get_one_info_int(const bcf_hdr_t *header, bcf1_t *record, const char *tag,int32_t &output);
    int bcf1_get_one_format_int(const bcf_hdr_t *header, bcf1_t *record, const char *tag,int32_t &output);

    //gets the index of a genotype likelihood for ploidy == 2
    int get_gl_index(int g0, int g1);

    //swaps the ath alle with the bth allele, rearranges PL/AD accordingly
    int bcf1_allele_swap(bcf_hdr_t *header, bcf1_t *record, int a,int b);

    //returns the string length of the right trimmed ref/alt (see https://academic.oup.com/bioinformatics/article/31/13/2202/196142)
    void right_trim(const char *ref,const char *alt,size_t &reflen,size_t &altlen);

    //Finds the index'th allele in from query in target. Retuns the allele index if found or -1 if missing.
    int find_allele(bcf1_t *target,bcf1_t *query,int index);

    //Adds the i'th allele from src to dst. FORMAT/INFO fields are unaffected.
    //Returns 0 if the allele was already present, 1 otherwise.
    int add_allele(bcf_hdr_t *hdr,bcf1_t *dst,bcf1_t *src,int index);

    //a heuristic for collapsing genotype likelihoods for variants that occur on two different rows but at the same position
    void collapse_gls(int ploidy,int num_alleles,std::vector< std::vector<int> > & pls,std::vector<int> & output);

    //this modifies the input vector, probably not the most efficient implementation of a median function
    float inplace_median(std::vector<int> & work); 

    //wrapper function for inplace_median. does not modify input
    float median(int *x, int n);

    //Fisher's exact test for per allele strand bias
    void fisher_sb_test(int *adf,int *adr,int num_allele,std::vector<float> & output,float maxret=1000.);

    std::string string_time();
    std::string generateUUID();


}
