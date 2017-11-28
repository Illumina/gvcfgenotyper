#ifndef GVCFMERGER_H
#define GVCFMERGER_H

#include <list>
#include <stdexcept>

extern "C" {
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/tbx.h>
}

#include "ggutils.hh"
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
    void genotype_homref_variant(int sample_index,DepthBlock & depth);
    void genotype_alt_variant(int sample_index,pair<std::deque<bcf1_t *>::iterator,std::deque<bcf1_t *>::iterator> & sample_variants);
    void genotype_sample(int sample_index);

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
    int32_t *_format_gt, *_format_gq, *_format_dp, *_format_dpf, *_format_ad, *_format_ps,*_format_pl,*_format_adf,*_format_adr,*_format_gqx;
    int32_t *_info_adf, *_info_adr, *_info_ac;
    int _num_pl,_mean_mq,_num_mq,_num_variants;
};

#endif

