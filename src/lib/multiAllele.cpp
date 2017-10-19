//
// Created by joconnell on 10/19/17.
//
#include <GVCFMerger.hpp>
#include <htslib/vcf.h>

inline bcf1_t *copy_alleles(bcf_hdr_t *hdr, bcf1_t *src)
{
    bcf_unpack(src,BCF_UN_ALL);
    assert(src->n_allele==2 || src->n_allele==3);//we only handle 1. binary alleles 2. binary alleles + symbolic 'X' allele
    if(src->n_allele==3)
    {
        assert(src->d.allele[2][0]=='X');
    }

    bcf1_t *alleles = bcf_init1();
    alleles->pos = src->pos;
    alleles->rid = src->rid;
    bcf_update_alleles(hdr,alleles,(const char**)src->d.allele,2);
    return(alleles);
}

multiAllele::multiAllele(int rid,int pos,bcf_hdr_t *hdr)
{
    _rid = rid;
    _pos = pos;
    _hdr = bcf_hdr_dup(hdr);
}

multiAllele::~multiAllele()
{
    for(auto rec=_records.begin();rec!=_records.end();rec++)
    {
        bcf_destroy(*rec);
    }
}

int multiAllele::allele(bcf1_t *record)
{
    if(record->rid!=_rid || record->pos!=_pos)
    {
        return(0);
    }

    bcf1_t *tmp = copy_alleles(record);
    auto location = std::find(_records.begin(),_records.end(),tmp);
    if(location==_records.end())
    {
        _records.push_back(tmp);
        return(_records.size());
    }
    else
    {
        bcf_destroy(tmp);
        return(location - _records.begin() + 1);
    }
}

bcf1_t *multiAllele::collapse()
{
    bcf1_t *ret = bcf_init1();
    
    return(ret);
}
