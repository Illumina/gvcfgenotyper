#ifndef GVCFMERGER_H
#define GVCFMERGER_H

#include "utils.hpp"
#include "GVCFReader.hpp"

extern "C" {
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/tbx.h>
}

class GVCFMerger
{
public:
    GVCFMerger(const vector<string> & input_files,const string &output_filename,const string & output_mode,const string & reference_genome,int buffer_size, const string& region="", const int is_file=0);
    ~GVCFMerger();
    void write_vcf();
    bcf1_t *next();
    
private:
    void build_header();
    bcf1_t *get_next_variant();
    void set_output_buffers_to_missing();
    vector<GVCFReader> _readers;
    int _num_gvcfs;
    bool all_readers_empty();
//stuff for output
    bcf1_t *_output_record;
    htsFile *_output_file;
    bcf_hdr_t *_output_header;
    int32_t *_format_gt,*_format_gq,*_format_dp,*_format_dpf,*_format_ad,*_format_ps,*format_pl;
};

#endif


