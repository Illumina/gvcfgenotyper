#include "kseq.h"
#include "kstring.h"
#include <deque>
#include <set>
#include <string>         // std::string
#include <locale>         // std::locale, std::toupper
#include <algorithm>
#include <iostream>
#include "needle.h"
#include "hts_utils.h"

KSTREAM_INIT(gzFile, gzread, 16384)


extern "C" {
#include "htslib/hts.h"
#include <htslib/vcf.h>
#include <htslib/tbx.h>
#include "vcfnorm.h"
}


#define ERR_REF_MISMATCH    -1
#define CHECK_REF_WARN 1

GVCFReader {

bool operator== (const bcf1_t *a,const bcf1_t *b)
{
    if(a->rid!=b->rid)
    {
	return(false);
    }
    else if(a->pos!=b->pos)
    {
	return(false);
    }
    else if(a->n_allele!=b->n_allele)
    {
	return(false);
    }
    else
    {
	for(int i=0;i<a->n_allele;i++)
	{
	    if(strcmp(a->d.allele[i],b->d.allele[i]))
	    {
		return(false);
	    }
	}
    }
    return(true);
}

bool operator!= (const bcf1_t *a,const bcf1_t *b)
{
    return(!a==b);
}

bool operator< (const bcf1_t *a,const bcf1_t *b)
{
    if(a->rid>b->rid)
    {
	return(false);
    }
    else if(a->pos>b->pos)
    {
	return(false);
    }
    else if(a->n_allele!=b->n_allele)
    {
	return(false);
    }
    else
    {
	for(int i=0;i<a->n_allele;i++)
	{
	    if(strcmp(a->d.allele[i],b->d.allele[i])>0)
	    {
		return(false);
	    }
	}
    }
    return(true);
}

bool operator> (const bcf1_t *a,const bcf1_t *b)
{
    return(!(a==b && a<b));
}


//small class to buffer bcf1_t records and sort them as they are inserted.
class VariantBuffer 
{
public:
    VarBuffer();
    ~VarBuffer();
    int push_back(bcf1_t *v);    //add a new variant (and sort if necessary)
    int flush_variant(int rid,int pos);//flush variants up to and including rid/pos
    int flush_buffer();//empty the buffer
    bool has_variant(bcf1_t *v);//does the buffer already have v?
  
private:
    int _last_pos,_ndup;
    deque<bcf1_t *> _buffer;  
    set <string> _seen; //list of seen variants at this position.
};

public:
    GVCFReader(const string & fname);
    int flush_buffer(int chrom,int pos);//empty buffer containing rows before and including chrom/pos
    bcf1_t *getCurrentRecord();
    
private:
    bcf_srs_t *_bcf_reader;
    bcf1_t *_current_bcf_record;
    int _num_duplicated_records;
}
