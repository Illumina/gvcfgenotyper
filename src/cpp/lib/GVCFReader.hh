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

#include "spdlog.h"


class GVCFReader
{
public:
    GVCFReader(const std::string &input_gvcf,Normaliser *normaliser, const int buffer_size,
               const string &region = "", const int is_file = 0);

    ~GVCFReader();

    int FlushBuffer();
    //empty buffer containing rows before and including chrom/pos
    int FlushBuffer(int chrom, int pos);
    //empty buffer containing rows before and including record
    int FlushBuffer(bcf1_t *record);

    pair<std::deque<bcf1_t *>::iterator,std::deque<bcf1_t *>::iterator> GetAllVariantsUpTo(bcf1_t *record);
    pair<std::deque<bcf1_t *>::iterator,std::deque<bcf1_t *>::iterator> GetAllVariantsInInterval(int chrom, int stop);
    bcf1_t *Front(); //return pointer to current vcf record
    bcf1_t *Pop(); //return pointer to current vcf record and remove it from buffer
    int ReadLines(const unsigned num_lines); //read at most num_lines
    size_t FillBuffer();

    //gets dp/dpf/gq (possibly interpolated) for a give interval
    void GetDepth(int rid, int start, int end, DepthBlock &db);
    bool IsEmpty();
    size_t GetNumVariants();
    size_t GetNumDepthBlocks();
    bcf_hdr_t *GetHeader();
    int ReadUntil(int rid, int pos);
    bool HasStrandAd();
    bool HasPl();
private:
    int _buffer_size;//ensure buffer has at least _buffer_size/2 variants avaiable (except at end of file)
    bcf_srs_t *_bcf_reader;//htslib synced reader.
    bcf1_t *_bcf_record;
    bcf_hdr_t *_bcf_header;
    VariantBuffer _variant_buffer;
    DepthBuffer _depth_buffer;
    Normaliser *_normaliser;
    std::shared_ptr<spdlog::logger> _lg;
};

#endif
