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

    bcf1_t *dst = bcf_init1();

    bcf_clear(dst);
    dst->rid  = src->rid;
    dst->pos  = src->pos;
    dst->rlen = src->rlen;
    dst->qual = src->qual;
    bcf_update_id(hdr,dst,".");
    kstring_t str = {0,0,nullptr};
    kputs(src->d.allele[0], &str);
    kputc(',', &str);
    kputs(src->d.allele[1],&str);
    bcf_update_alleles_str(hdr,dst,str.s);
    free(str.s);
    return(dst);
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

    bcf1_t *tmp = copy_alleles(_hdr,record);
    auto location = _records.begin();
    while(location!=_records.end() && bcf1_not_equal(*location,tmp))
    {
        location++;
    }
    if(location==_records.end())
    {
        _records.push_back(tmp);
        return((int)_records.size());
    }
    else
    {
        bcf_destroy(tmp);
        return((int) std::distance(_records.begin(),location)+1);
    }
}

bcf1_t *multiAllele::get_max()
{
    bcf1_t *ret = nullptr;
    for(auto rec=_records.begin();rec!=_records.end();rec++)
    {
        if(ret== nullptr || bcf1_greater_than(*rec,ret))
        {
//            print_variant(*rec);
            ret = *rec;
        }
    }
    return(ret);
}

void multiAllele::collapse(bcf1_t *output)
{
    assert(!_records.empty());
    output->pos = _pos;
    output->rid = _rid;

    int num_alleles = (int)_records.size()+1;
    char **new_alleles = (char **)malloc(sizeof(char*) * num_alleles);
    char *ptr_to_longest= nullptr;
    size_t max_ref_len = 0;
    for(auto rec=_records.begin();rec!=_records.end();rec++)
    {
        size_t curlen = strlen((*rec)->d.allele[0]);
        if(curlen>max_ref_len)
        {
            ptr_to_longest = (*rec)->d.allele[0];
            max_ref_len=curlen;
        }
    }
    assert(ptr_to_longest!=nullptr);
    int index=0;
    new_alleles[index] = (char *)malloc(max_ref_len+1);
    strcpy(new_alleles[index++],ptr_to_longest);
    for(auto rec=_records.begin();rec!=_records.end();rec++)
    {
        if(strncmp(new_alleles[0],(*rec)->d.allele[0],strlen((*rec)->d.allele[0]))!=0)
        {
            std::cerr << _pos+1 << " " <<  new_alleles[0] << "!=" <<(*rec)->d.allele[0] << std::endl;
            print_variant(*rec);
            throw std::runtime_error("inconsistent REF on first allele");
        }
        size_t old_ref_len = strlen((*rec)->d.allele[0]);
        size_t rightpad = max_ref_len - old_ref_len;
        char *new_alt=(*rec)->d.allele[1];
        new_alleles[index]=(char *)malloc(strlen(new_alt)+rightpad+1);
        memcpy(new_alleles[index],new_alt,strlen(new_alt));
        memcpy(new_alleles[index++]+strlen(new_alt),new_alleles[0]+old_ref_len,rightpad+1);
    }
    bcf_update_alleles(_hdr,output,(const char**)new_alleles,num_alleles);
    for(size_t i=0;i<_records.size();i++)
    {
        free(new_alleles[i]);
    }
    free(new_alleles);
}
