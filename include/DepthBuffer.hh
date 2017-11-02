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

    void push_back(DepthBlock db);

    DepthBlock *pop();

    DepthBlock *back();

    DepthBlock *front();

    DepthBlock intersect(const DepthBlock &db);

    int flush_buffer();

    int flush_buffer(int rid, int pos);

    int interpolate(int rid, int start, int end, DepthBlock &db);//interpolates depth for an interval a<=x<
    size_t size();

private:
    deque<DepthBlock> _buffer;
};

#endif //GVCFGENOTYPER_DEPTHBUFFER_HH
