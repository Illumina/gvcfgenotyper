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

    void Init(bcf_hdr_t *hdr);
    void SetPosition(int rid, int pos);
    int Allele(bcf1_t *record,int index=1);
    int AlleleIndex(bcf1_t *record,int index);
    void Collapse(bcf1_t *output);
    int GetNumAlleles() {return _records.size();};
    bcf1_t *GetMax();//returns the maximum allele (as defined by bcf1_t_less_than)
    int Clear();//wipes the _records

    //accessors/mutators
    int pos() {return _pos;};
    int rid() {return _rid;};
    void print();

private:
    int _rid,_pos;
    bcf_hdr_t *_hdr;
    list<bcf1_t*> _records;//FIXME. we should probably use some sorted data structure + binary search here. in practice in might not matter.
};


#endif //GVCFGENOTYPER_MULTIALLELE_HH
