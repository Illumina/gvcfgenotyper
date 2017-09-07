#ifndef GVCFREADER_H
#define GVCFREADER_H

#include <deque>
#include <set>
#include <string>         // std::string
#include <locale>         // std::locale, std::toupper
#include <algorithm>
#include <iostream>
#include "utils.hpp"
//#include "hts_utils.h"

extern "C" {
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/tbx.h>
#include "vcfnorm.h"
}

#define ERR_REF_MISMATCH    -1
#define CHECK_REF_WARN 1


#define MROWS_SPLIT 1
#define MROWS_MERGE  2
//this basically wraps bcftools norm in a class.
//TODO: investigate replacing this with invariant components
class Normaliser {
public:
    Normaliser(const std::string & ref_fname,bcf_hdr_t *hdr);
    ~Normaliser();
    std::vector<bcf1_t *> atomise(bcf1_t *rec);
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
    int flush_buffer(int rid,int pos);//flush variants up to and including rid/pos
    int flush_buffer();//empty the buffer
    bool has_variant(bcf1_t *v);//does the buffer already have v?
    bcf1_t *front(); //return pointer to current vcf record
    bcf1_t *pop(); //return pointer to current vcf record and remove it from buffer
    bool empty();  
    size_t size();
private:
    int _last_pos, _num_duplicated_records;
    deque<bcf1_t *> _buffer;  
    set <std::string> _seen; //list of seen variants at this position.
};

class GVCFReader {
public:
    GVCFReader(const std::string & input_gvcf,const std::string & reference_genome_fasta,int buffer_size);
    ~GVCFReader();
    int flush_buffer();
    int flush_buffer(int chrom,int pos);//empty buffer containing rows before and including chrom/pos
    bcf1_t *front(); //return pointer to current vcf record
    bcf1_t *pop(); //return pointer to current vcf record and remove it from buffer
    int read_lines(int num_lines); //read at most num_lines
    int fill_buffer(int num_lines);
    bool empty();
    const bcf_hdr_t *getHeader();
private:
    int  _buffer_size;//ensure buffer has at least _buffer_size/2 variants avaiable (except at end of file)
    bcf_srs_t *_bcf_reader;//htslib synced reader.
    bcf1_t *_bcf_record;
    bcf_hdr_t *_bcf_header;
    VariantBuffer _variant_buffer;
    Normaliser *_normaliser;        
};

#endif