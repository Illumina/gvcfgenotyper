//
// Created by O'Connell, Jared on 11/2/17.
//

#ifndef GVCFGENOTYPER_DEPTHBUFFER_HH
#define GVCFGENOTYPER_DEPTHBUFFER_HH

extern "C" {
#include <htslib/vcf.h>
}

#include "DepthBlock.hh"

class DepthBuffer
{
public:
    DepthBuffer()
    {};

    ~DepthBuffer()
    {};

    // performs additional check on db, if passed moves db into internal buffer
    // which invalidates db
    void push_back(const DepthBlock& db);

    DepthBlock *Pop();
    DepthBlock *Back();
    //DepthBlock *Front();
    DepthBlock Intersect(const DepthBlock &db);
    int FlushBuffer();
    int FlushBuffer(const int rid, const int pos);
    int Interpolate(const int rid, const int start, const int end, DepthBlock &db);//interpolates depth for an interval a<=x<

    //accessors/mutators
    size_t Size();

    bool Empty();

    void dump();

private:
    deque<DepthBlock> _buffer;
};

#endif //GVCFGENOTYPER_DEPTHBUFFER_HH
