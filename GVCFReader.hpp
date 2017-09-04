#include <deque>
#include <set>
#include <string>         // std::string
#include <locale>         // std::locale, std::toupper
#include <algorithm>
#include <iostream>
#include "utils.h"
#include "BCFHelpers.hpp"
//#include "hts_utils.h"


extern "C" {
#include "htslib/hts.h"
#include <htslib/vcf.h>
#include <htslib/tbx.h>
#include "vcfnorm.h"
}


#define ERR_REF_MISMATCH    -1
#define CHECK_REF_WARN 1


bool operator== (const bcf1_t & a,const bcf1_t & b)
{
    if(a.rid!=b.rid)
    {
	return(false);
    }
    else if(a.pos!=b.pos)
    {
	return(false);
    }
    else if(a.n_allele!=b.n_allele)
    {
	return(false);
    }
    else
    {
	for(int i=0;i<a.n_allele;i++)
	{
	    if(strcmp(a.d.allele[i],b.d.allele[i]))
	    {
		return(false);
	    }
	}
    }
    return(true);
}

bool operator!= (const bcf1_t & a,const bcf1_t & b)
{
    return(!(a==b));
}

bool operator< (const bcf1_t & a,const bcf1_t & b)
{
    if(a.rid>b.rid)
    {
	return(false);
    }
    else if(a.pos>b.pos)
    {
	return(false);
    }
    else if(a.n_allele!=b.n_allele)
    {
	return(false);
    }
    else
    {
	for(int i=0;i<a.n_allele;i++)
	{
	    if(strcmp(a.d.allele[i],b.d.allele[i])>0)
	    {
		return(false);
	    }
	}
    }
    return(true);
}

bool operator> (const bcf1_t & a,const bcf1_t & b)
{
    return(!(a==b && a<b));
}



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
  
private:
    int _last_pos, _num_duplicated_records;    
    deque<bcf1_t *> _buffer;  
    set <std::string> _seen; //list of seen variants at this position.
};

#define MROWS_SPLIT 1
#define MROWS_MERGE  2
//this basically wraps bcftools norm in a class.
//TODO: investigate replacing this with invariant components
class Normaliser {
public:
    Normaliser(const std::string & ref_fname);
    ~Normaliser();
    std::vector<bcf1_t *> atomise(bcf1_t *rec,bcf_hdr_t *hdr);
private:
    args_t *_norm_args;
};


class GVCFReader {
public:
    GVCFReader(const std::string & input_gvcf,const std::string & reference_genome_fasta);
    ~GVCFReader();
    int flush_buffer(int chrom,int pos);//empty buffer containing rows before and including chrom/pos
    bcf1_t *get_current_record(); //return pointer to current vcf record
    int read_lines(int num_lines); //read at most num_lines

private:
    bcf_srs_t *_bcf_reader;//htslib synced reader.
    bcf1_t *_bcf_record;
    bcf_hdr_t *_bcf_header;
    VariantBuffer _variant_buffer;
    Normaliser *_normaliser;        
};
