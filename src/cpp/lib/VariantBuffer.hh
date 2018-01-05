//
// Created by O'Connell, Jared on 11/2/17.
//

#ifndef GVCFGENOTYPER_VARIANTBUFFER_HH
#define GVCFGENOTYPER_VARIANTBUFFER_HH

#include <deque>

extern "C" {
#include <htslib/vcf.h>
}

#include "ggutils.hh"

//small class to buffer bcf1_t records and sort them as they are inserted.
class VariantBuffer
{
public:
    VariantBuffer();
    ~VariantBuffer();

    //add a new variant (and sort if necessary), needs header for duplicate check
    // Warning: if the variant already occurs in the buffer, bcf_destroy is called on it
    int PushBack(const bcf_hdr_t *hdr, bcf1_t *v);
    int FlushBuffer(int rid, int pos);//flush variants up to and including rid/pos
    int FlushBuffer();//empty the buffer
    int FlushBuffer(bcf1_t *record);
    pair<std::deque<bcf1_t *>::iterator,std::deque<bcf1_t *>::iterator> GetAllVariantsInInterval(int chrom, int stop);//gets all variants in interval start<=x<=stop
    pair<std::deque<bcf1_t *>::iterator,std::deque<bcf1_t *>::iterator> GetAllVariantsUpTo(bcf1_t *record);//gets all variants in interval start<=x<=stop

    bool HasVariant(const bcf_hdr_t *hdr, bcf1_t *v);//does the buffer already have v? Swap info fields if duplicate is hom-ref
    bcf1_t *Front(); //return pointer to current vcf record
    bcf1_t *Back(); //return pointer to last vcf record
    bcf1_t *Pop(); //return pointer to current vcf record and remove it from buffer
    bool IsEmpty();

    size_t Size();
    size_t GetNumDuplicatedRecords() const { return _num_duplicated_records;};

private:
    size_t  _num_duplicated_records;
    deque<bcf1_t *> _buffer;
    set<std::string> _seen; //list of seen variants at this position.
};

#endif //GVCFGENOTYPER_VARIANTBUFFER_HH
