#include "GVCFReader.hpp"

DepthBlock::DepthBlock()
{
    set_missing();
}

DepthBlock::DepthBlock(int rid,int start,int end,int dp,int dpf,int gq)
{
    assert(end>=0 && start>=0);
    assert(end>=start);
    _rid=rid;
    _start=start;
    _end=end;
    _dp=dp;
    _dpf=_dpf;
    _gq=gq;
}

void DepthBlock::set_missing()
{
    _dp = _gq = _dpf = bcf_int32_missing;
}

void DepthBlock::zero()
{
    _dp = _gq = _dpf = 0;
}

void DepthBlock::add(const DepthBlock & db)
{
    _dp += db._dp;
    _dpf += db._dpf;
    _gq += db._gq;
}

void DepthBlock::divide(int n)
{
    if(n>1)
    {
	_dp = (int)round((float)_dp/(float)n);
	_gq =  (int)round((float)_gq/(float)n);
	_dpf =  (int)round((float)_dpf/(float)n);
    }
}

int DepthBlock::intersect_size(int rid,int start,int stop)
{
    if(_rid!=rid)
    {
	return(0);
    }
    if(_end < start || _start>stop)
    {
	return(0);
    }
    return(min(stop,_end) - max(_start,start) + 1);    
}

int DepthBlock::intersect_size(const DepthBlock &db)
{
    return(intersect_size(db._rid,db._start,db._end) );
}

