//
// Created by O'Connell, Jared on 11/1/17.
//

#ifndef GVCFGENOTYPER_MULTIALLELE_HH
#define GVCFGENOTYPER_MULTIALLELE_HH


#include <list>
#include <stdexcept>

extern "C" {
#include <htslib/vcf.h>
}

#include "ggutils.hh"

//gathers multiple alleles. kind of like a set() for bcf1_t
class multiAllele
{
public:
    multiAllele();
    ~multiAllele();

    void init(bcf_hdr_t *hdr);
    void setPosition(int rid,int pos);
    int allele(bcf1_t *record);
    void collapse(bcf1_t *output);
    int get_pos() {return _pos;};
    int get_rid() {return _rid;};
    int num_alleles() {return _records.size();};
    bcf1_t *get_max();//returns the maximum allele (as defined by bcf1_t_less_than)
    int clear();//wipes the _records
private:
    int _rid,_pos;
    bcf_hdr_t *_hdr;
    list<bcf1_t*> _records;//FIXME. we should probably use some sorted data structure + binary search here. in practice in might not matter.
};


#endif //GVCFGENOTYPER_MULTIALLELE_HH
