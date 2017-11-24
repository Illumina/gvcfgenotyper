#pragma once

#define __STDC_LIMIT_MACROS

#include <stdint.h>
#include "math.h"
#include <stdlib.h>

#include <limits>
#include <set>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>

using namespace std;

extern "C" {
#include <htslib/bgzf.h>
#include <htslib/tbx.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <htslib/synced_bcf_reader.h>
}


namespace ggutils
{
    int read_text_file(const string &fname, vector<string> &output);

    void print_variant(bcf_hdr_t const *header, bcf1_t *record);

    void print_variant(bcf1_t *record);

    int *zeros(int n);

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

    bool is_complex(bcf1_t *record);

    int get_variant_rank(bcf1_t *record);

    int get_end_of_gvcf_block(bcf_hdr_t *header, bcf1_t *record);

    int get_ploidy(bcf_hdr_t *header, bcf1_t *record);

    int get_end_of_variant(bcf1_t *record);

    bool bcf1_equal(bcf1_t *a, bcf1_t *b);

    bool bcf1_less_than(bcf1_t *a, bcf1_t *b);

    bool bcf1_greater_than(bcf1_t *a, bcf1_t *b);

    bool bcf1_leq(bcf1_t *a, bcf1_t *b);

    bool bcf1_geq(bcf1_t *a, bcf1_t *b);

    bool bcf1_not_equal(bcf1_t *a, bcf1_t *b);

    size_t get_number_of_likelihoods(int ploidy, int num_allele);

    int phred(float l);

    float unphred(int pl);

    int factorial(int x);

    int choose(int n, int k);

    //these are simple htslib wrappers shortcuts to get scalar ints/floats (not appropriate for arrays)
    float bcf1_get_one_info_float(bcf_hdr_t *header, bcf1_t *record, const char *tag);
    float bcf1_get_one_format_float(bcf_hdr_t *header, bcf1_t *record, const char *tag);
    int bcf1_get_one_info_int(bcf_hdr_t *header, bcf1_t *record, const char *tag);
    int bcf1_get_one_format_int(bcf_hdr_t *header, bcf1_t *record, const char *tag);

    //gets the index of a genotype likelihood for ploidy == 2
    int get_gl_index(int g0, int g1);
}