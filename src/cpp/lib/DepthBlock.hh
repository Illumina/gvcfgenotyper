//
// Created by O'Connell, Jared on 11/2/17.
//

#ifndef GVCFGENOTYPER_DEPTHBLOCK_HH
#define GVCFGENOTYPER_DEPTHBLOCK_HH

extern "C" {
#include <htslib/vcf.h>
}

#include "ggutils.hh"

//simple class that stores the pertinent values from a GVCF homref block
class DepthBlock
{
public:
    DepthBlock();

    DepthBlock(int rid, int start, int end, int dp, int dpf, int gq,int ploidy);

    inline bool operator == (const DepthBlock& db) const {
        return (_rid==db._rid &&
                _start==db._start &&
                _end==db._end &&
                _dp==db._dp &&
                _dpf==db._dpf &&
                _gq==db._gq
               );
    }

    DepthBlock Intersect(const DepthBlock &db);
    DepthBlock Intersect(int rid, int start, int end);
    int IntersectSize(int rid, int a, int b) const;
    int IntersectSize(const DepthBlock &db) const;

    void SetToMissing();//set all values to bcftools missing
    void SetToZero();//zero all values
    void Add(const DepthBlock &db);
    void Divide(int n);

    int rid() const { return _rid; }
    int start() const { return _start; }
    int end() const { return _end; }
    int dp() const { return _dp; }
    int gq() const { return _gq; }
    int dpf() const { return _dpf; }
    int ploidy();
    int size() const;

private:
    int _rid, _start, _end, _dp, _dpf, _gq, _ploidy;
};

#endif //GVCFGENOTYPER_DEPTHBLOCK_HH


