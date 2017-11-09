#ifndef GVCFMERGER_H
#define GVCFMERGER_H

#include <list>
#include <stdexcept>

extern "C" {
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/tbx.h>
}

#include "utils.hh"
#include "GVCFReader.hh"
#include "multiAllele.hh"
#include "Genotype.hh"

class GVCFMerger
{
public:
    GVCFMerger(const vector<string> &input_files, const string &output_filename, const string &output_mode,
               const string &reference_genome, int buffer_size, const string &region = "", const int is_file = 0);
    ~GVCFMerger();
    void write_vcf();
    bcf1_t *next();
    int get_next_variant();

private:
    void build_header();
    void set_output_buffers_to_missing(int num_alleles);
    vector<GVCFReader> _readers;
    size_t _num_gvcfs;
    bool all_readers_empty();
    multiAllele _record_collapser;

//stuff for output
    bcf1_t *_output_record;
    htsFile *_output_file;
    bcf_hdr_t *_output_header;
    int32_t *_format_gt, *_format_gq, *_format_dp, *_format_dpf, *_format_ad, *_format_ps,*_format_pl;
    int _num_pl;
    int _num_variants;
};

#endif

