//
// Created by O'Connell, Jared on 11/2/17.
//

#ifndef GVCFGENOTYPER_VARIANTBUFFER_HH
#define GVCFGENOTYPER_VARIANTBUFFER_HH

#include <deque>

extern "C" {
#include <htslib/vcf.h>
}

#include "utils.hpp"

//small class to buffer bcf1_t records and sort them as they are inserted.
class VariantBuffer
{
public:
    VariantBuffer();
    ~VariantBuffer();

    int push_back(bcf1_t *v);    //add a new variant (and sort if necessary)
    int flush_buffer(int rid, int pos);//flush variants up to and including rid/pos
    int flush_buffer();//empty the buffer
    int flush_buffer(bcf1_t *record);
    pair<std::deque<bcf1_t *>::iterator,std::deque<bcf1_t *>::iterator> get_all_variants_in_interval(int chrom,int stop);//gets all variants in interval start<=x<=stop
    pair<std::deque<bcf1_t *>::iterator,std::deque<bcf1_t *>::iterator> get_all_variants_up_to(bcf1_t *record);//gets all variants in interval start<=x<=stop

    bool has_variant(bcf1_t *v);//does the buffer already have v?
    bcf1_t *front(); //return pointer to current vcf record
    bcf1_t *back(); //return pointer to last vcf record
    bcf1_t *pop(); //return pointer to current vcf record and remove it from buffer
    bool empty();

    size_t size();

private:
    int  _num_duplicated_records;
    deque<bcf1_t *> _buffer;
    set<std::string> _seen; //list of seen variants at this position.
};

#endif //GVCFGENOTYPER_VARIANTBUFFER_HH
