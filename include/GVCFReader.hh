#ifndef GVCFREADER_H
#define GVCFREADER_H

#include <set>
#include <string>         // std::string
#include <locale>         // std::locale, std::toupper
#include <algorithm>
#include <iostream>
#include "numeric"

extern "C" {
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include "vcfnorm.h"
}

#include "ggutils.hh"
#include "Normaliser.hh"
#include "VariantBuffer.hh"
#include "DepthBuffer.hh"


class GVCFReader
{
public:
    GVCFReader(const std::string &input_gvcf, const std::string &reference_genome_fasta, const int buffer_size,
               const string &region = "", const int is_file = 0);

    ~GVCFReader();

    int flush_buffer();

    int flush_buffer(int chrom, int pos);//empty buffer containing rows before and including chrom/pos
    int flush_buffer(bcf1_t *record);

    pair<std::deque<bcf1_t *>::iterator,std::deque<bcf1_t *>::iterator> get_all_variants_up_to(bcf1_t *record);
    pair<std::deque<bcf1_t *>::iterator,std::deque<bcf1_t *>::iterator> get_all_variants_in_interval(int chrom,int stop);
    bcf1_t *front(); //return pointer to current vcf record
    bcf1_t *pop(); //return pointer to current vcf record and remove it from buffer
    int read_lines(const unsigned num_lines); //read at most num_lines
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
