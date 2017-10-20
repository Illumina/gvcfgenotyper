#ifndef GVCFMERGER_H
#define GVCFMERGER_H

#include "utils.hpp"
#include "GVCFReader.hpp"
#include <list>


extern "C" {
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/tbx.h>
}


//gathers multiple alleles. kind of like a set() for bcf1_t
class multiAllele
{
public:
    multiAllele(int rid,int pos,bcf_hdr_t *hdr);
    ~multiAllele();
    int allele(bcf1_t *record);
    void collapse(bcf1_t *output);
    int get_pos() {return _pos;};
    int get_rid() {return _rid;};

private:
    int _rid,_pos;
    bcf_hdr_t *_hdr;
    list<bcf1_t*> _records;//FIXME. we should probably use some sorted data structure + binary search here. in practice in might not matter.
};


class GVCFMerger
{
public:
    GVCFMerger(const vector<string> &input_files, const string &output_filename, const string &output_mode,
               const string &reference_genome, int buffer_size, const string &region = "", const int is_file = 0);
    ~GVCFMerger();
    void write_vcf();
    bcf1_t *next();
    multiAllele get_next_variant();

private:
    void build_header();
    void set_output_buffers_to_missing();
    vector<GVCFReader> _readers;
    int _num_gvcfs;
    bool all_readers_empty();

//stuff for output
    bcf1_t *_output_record;
    htsFile *_output_file;
    bcf_hdr_t *_output_header;
    int32_t *_format_gt, *_format_gq, *_format_dp, *_format_dpf, *_format_ad, *_format_ps, *format_pl;
    bool var_without_gq_seen;
};

#endif

