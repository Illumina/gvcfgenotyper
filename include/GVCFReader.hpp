#ifndef GVCFREADER_H
#define GVCFREADER_H

#include <deque>
#include <set>
#include <string>         // std::string
#include <locale>         // std::locale, std::toupper
#include <algorithm>
#include <iostream>
#include "utils.hpp"
#include "numeric"

extern "C" {
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/tbx.h>
#include "vcfnorm.h"
}

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

class Genotype
{
public:
    Genotype(bcf_hdr_t *header, bcf1_t *record);

    Genotype(int ploidy, int num_allele);

    Genotype marginalise(int index);

    ~Genotype();

    void setDepthFromAD();

    int update_bcf1_t(bcf_hdr_t *header, bcf1_t *record);

    int *_gt, *_ad, *_gq, *_dp, *_dpf, *_pl;
    int _num_allele, _num_pl, _ploidy, _num_gt, _num_ad, _num_gq, _num_dp, _num_dpf, _num_gl;
    std::vector<float> _gl;
};

//this basically wraps bcftools norm in a class.
//TODO: investigate replacing this with invariant components
class Normaliser
{
public:
    Normaliser(const std::string &ref_fname, bcf_hdr_t *hdr);

    ~Normaliser();

    std::vector<bcf1_t *> unarise(bcf1_t *rec);

private:
    args_t *_norm_args;
    bcf_hdr_t *_hdr;
};

//small class to buffer bcf1_t records and sort them as they are inserted.
class VariantBuffer
{
public:
    VariantBuffer();

    ~VariantBuffer();

    int push_back(bcf1_t *v);    //add a new variant (and sort if necessary)
    int flush_buffer(int rid, int pos);//flush variants up to and including rid/pos
    int flush_buffer();//empty the buffer
    int flush_buffer(const bcf1_t *record);
    vector<bcf1_t *> get_all_variants_in_interval(int chrom,int stop);//gets all variants in interval start<=x<=stop

    bool has_variant(bcf1_t *v);//does the buffer already have v?
    bcf1_t *front(); //return pointer to current vcf record
    bcf1_t *back(); //return pointer to last vcf record
    bcf1_t *pop(); //return pointer to current vcf record and remove it from buffer
    bool empty();

    size_t size();

private:
    int _last_pos, _num_duplicated_records;
    deque<bcf1_t *> _buffer;
    set<std::string> _seen; //list of seen variants at this position.
};


//simple class that stores the pertinent values from a GVCF homref block
class DepthBlock
{
public:
    DepthBlock();

    DepthBlock(int rid, int start, int end, int dp, int dpf, int gq);

    DepthBlock intersect(const DepthBlock &db);

    DepthBlock intersect(int rid, int start, int end);

    int intersect_size(int rid, int a, int b) const;

    int intersect_size(const DepthBlock &db) const;

    int size() const;

    void set_missing();//set all values to bcftools missing
    void zero();//zero all values
    void add(const DepthBlock &db);

    void divide(int n);

    int _rid, _start, _end, _dp, _gq, _dpf;
};

class DepthBuffer
{
public:
    DepthBuffer()
    {};

    ~DepthBuffer()
    {};

    void push_back(DepthBlock db);

    DepthBlock *pop();

    DepthBlock *back();

    DepthBlock *front();

    DepthBlock intersect(const DepthBlock &db);

    int flush_buffer();

    int flush_buffer(int rid, int pos);

    int interpolate(int rid, int start, int end, DepthBlock &db);//interpolates depth for an interval a<=x<
    size_t size();

private:
    deque<DepthBlock> _buffer;
};

class GVCFReader
{
public:
    GVCFReader(const std::string &input_gvcf, const std::string &reference_genome_fasta, const int buffer_size,
               const string &region = "", const int is_file = 0);

    ~GVCFReader();

    int flush_buffer();

    int flush_buffer(int chrom, int pos);//empty buffer containing rows before and including chrom/pos
    int flush_buffer(const bcf1_t *record);

    vector<bcf1_t *> get_all_variants_in_interval(int chrom,int stop);
    bcf1_t *front(); //return pointer to current vcf record
    bcf1_t *pop(); //return pointer to current vcf record and remove it from buffer
    int read_lines(int num_lines); //read at most num_lines
    size_t fill_buffer();

    void
    get_depth(int rid, int start, int end, DepthBlock &db);//gets dp/dpf/gq (possibly interpolated) for a give interval
    bool empty();

    size_t get_num_variants();

    size_t get_num_depth();

    const bcf_hdr_t *get_header();

    int read_until(int rid, int pos);

private:
    int _buffer_size;//ensure buffer has at least _buffer_size/2 variants avaiable (except at end of file)
    bcf_srs_t *_bcf_reader;//htslib synced reader.
    bcf1_t *_bcf_record;
    bcf_hdr_t *_bcf_header;
    VariantBuffer _variant_buffer;
    DepthBuffer _depth_buffer;
    Normaliser *_normaliser;
};

#endif
