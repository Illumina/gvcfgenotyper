//
// Created by joconnell on 10/19/17.
//
#include <GVCFMerger.hh>
#include <htslib/vcf.h>

inline bcf1_t *copy_alleles(bcf_hdr_t *hdr, bcf1_t *src,int index)
{
    assert(src->n_allele>1);
    if(!(index>0 && index<src->n_allele)) {
        ggutils::print_variant(hdr,src);
        cerr << "bad index = " << index << "\n";
        //ggutils::die("bad index = "+std::to_string(index) + " nal=" + std::to_string(src->n_allele));
        assert(0);
    }
    bcf_unpack(src,BCF_UN_ALL);
    size_t rlen,alen;
    ggutils::right_trim(src->d.allele[0],src->d.allele[index], rlen,alen);

    bcf1_t *dst = bcf_init();
    bcf_clear(dst);
    dst->rid  = src->rid;
    dst->pos  = src->pos;

    kstring_t str = {0,0,nullptr};
    kputsn(src->d.allele[0],rlen,&str);
    kputc(',', &str);
    kputsn(src->d.allele[index],alen,&str);
    bcf_update_alleles_str(hdr,dst,str.s);
    free(str.s);
    return(dst);
}

multiAllele::multiAllele()
{
    _rid = -1;
    _pos = -1;
    _hdr = nullptr;
}

void multiAllele::SetPosition(int rid, int pos)
{
    assert(_hdr!= nullptr);
    Clear();
    _rid = rid;
    _pos = pos;
}

void multiAllele::Init(bcf_hdr_t *hdr)
{
    _hdr = bcf_hdr_dup(hdr);
}

int multiAllele::Clear()
{
    _rid = -1;
    _pos = -1;
    for (auto rec = _records.begin(); rec != _records.end(); rec++)
    {
        bcf_destroy(*rec);
    }
    int ret = _records.size();
    _records.clear();
    return(ret);
}

multiAllele::~multiAllele()
{
    Clear();
    bcf_hdr_destroy(_hdr);
}

int multiAllele::Allele(bcf1_t *record,int index)
{
    assert(_hdr!=nullptr);
    assert(_rid>=0 && _pos>=0);
    if(record->rid!=_rid || record->pos!=_pos) return(0);

    bcf1_t *tmp = copy_alleles(_hdr,record,index);
    auto location = _records.begin();
    while(location!=_records.end() && ggutils::bcf1_not_equal(*location,tmp))
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

int multiAllele::AlleleIndex(bcf1_t *record,int index)
{
    if(index==0) return 0;
    assert(_hdr!=nullptr);
    assert(_rid>=0 && _pos>=0);
    if(record->rid!=_rid || record->pos!=_pos)
    {
        return(0);
    }

    bcf1_t *tmp = copy_alleles(_hdr,record,index);
    auto location = _records.begin();
    //while(location!=_records.end() && ggutils::bcf1_not_equal(*location,tmp)) location++;
    while(location!=_records.end() && ggutils::find_allele(*location,record,index)!=1) location++;

    if(location==_records.end())
    {
        print();
        ggutils::print_variant(record);
        ggutils::die("multiAllele: variant "+std::to_string(index)+" not found");
        return(-1);
    }
    else
    {
        bcf_destroy(tmp);
        return((int) std::distance(_records.begin(),location)+1);
    }
}

bcf1_t *multiAllele::GetMax()
{
    bcf1_t *ret = nullptr;
    for(auto rec=_records.begin();rec!=_records.end();rec++)    {

        if(ret== nullptr ||  ggutils::bcf1_greater_than(*rec,ret))
        {
            //ggutils::print_variant(*rec);
            ret = *rec;
        }
    }
    return(ret);
}

void multiAllele::Collapse(bcf1_t *output)
{
    assert(!_records.empty());
    output->pos = _pos;
    output->rid = _rid;

    size_t num_alleles = _records.size()+1;
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
            ggutils::print_variant(*rec);
            ggutils::die("inconsistent REF on first allele");
        }
        size_t old_ref_len = strlen((*rec)->d.allele[0]);
        size_t rightpad = max_ref_len - old_ref_len;
        char *new_alt=(*rec)->d.allele[1];
        new_alleles[index]=(char *)malloc(strlen(new_alt)+rightpad+1);
        memcpy(new_alleles[index],new_alt,strlen(new_alt));
        memcpy(new_alleles[index++]+strlen(new_alt),new_alleles[0]+old_ref_len,rightpad+1);
    }
    bcf_update_alleles(_hdr,output,(const char**)new_alleles,num_alleles);
    for(size_t i=0;i<num_alleles;i++) free(new_alleles[i]);
    free(new_alleles);
}

void multiAllele::print()
{
    std::cerr<<"multiAllele start"<<std::endl;
    for(auto it=_records.begin();it!=_records.end();it++)
    {
        ggutils::print_variant(_hdr,*it);
        std::cerr<<"rank="<<ggutils::get_variant_rank(*it)<<std::endl;
    }

    std::cerr<<"multiAllele end"<<std::endl;
}
