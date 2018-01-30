
#ifndef GVCFMERGER_H
#define GVCFMERGER_H

#include <list>
#include <stdexcept>

#include "spdlog.h"


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
    GVCFMerger(const vector<string> &input_files,
	           const string &output_filename,
	           const string &output_mode,
               const string &reference_genome,
	           int buffer_size,
	           const string &region = "",
               const int is_file = 0,
               bool ignore_non_matching_ref=false,
               bool force_samples=false);
    ~GVCFMerger();
    void write_vcf();
    bcf1_t *next();
    int GetNextVariant();
    void SetMaxAlleles(size_t max_alleles) {_max_alleles=max_alleles;};

private:
    void GenotypeHomrefVariant(int sample_index, DepthBlock &depth);
    void GenotypeAltVariant(int sample_index,bcf1_t *sample_variants);
    void GenotypeSample(int sample_index);
    void UpdateFormatAndInfo();
    void BuildHeader();
    void SetOutputBuffersToMissing(int num_alleles);
    bool AreAllReadersEmpty();


    multiAllele _record_collapser;
    vector<GVCFReader> _readers;
    size_t _num_gvcfs;
    bcf1_t *_output_record;
    htsFile *_output_file;
    bcf_hdr_t *_output_header;
    ggutils::vcf_data_t *format;//stores all our format fields.
    int32_t *_info_adf, *_info_adr, *_info_ac;
    int _mean_mq,_num_mq,_num_variants;
    size_t _num_ps_written;
    bool _has_strand_ad,_has_pl;
    Normaliser *_normaliser;
    std::shared_ptr<spdlog::logger> _lg;
    bool _force_samples;
	size_t _max_alleles;
};

#endif

