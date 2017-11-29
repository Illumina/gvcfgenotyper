#include "DepthBlock.hh"

#include<iostream>

DepthBlock::DepthBlock()
    : _rid(0),_start(0),_end(0)
{
    set_missing();
}

DepthBlock::DepthBlock(int rid, int start, int end, int dp, int dpf, int gq, int ploidy)
{
    assert(end >= 0 && start >= 0);
    assert(end >= start);
    _rid = rid;
    _start = start;
    _end = end;
    _dp = dp;
    _dpf = dpf;
    _gq = gq;
    _ploidy = ploidy;
}

void DepthBlock::set_missing()
{
    _dp = _gq = _dpf = bcf_int32_missing;
}

void DepthBlock::zero()
{
    _dp = _gq = _dpf = 0;
}

void DepthBlock::add(const DepthBlock &db)
{
    if (_dp == bcf_int32_missing || _dpf == bcf_int32_missing || _gq == bcf_int32_missing)
    {
        _dp = db._dp;
        _dpf = db._dpf;
        _gq = db._gq;
        _rid = db._rid;
        _start = db._start;
        _end = db._end;
    }
    else
    {
        assert(db._rid == _rid && db._start == (_end + 1));
        float length1 = size();
        float length2 = db.size();
        float p1 = length1 / (length1 + length2);
        float p2 = length2 / (length1 + length2);
        _end = db._end;
        _dp = round(p1 * _dp + p2 * db._dp);
        _dpf = round(p1 * _dpf + p2 * db._dpf);
        _gq = round(p1 * _gq + p2 * db._gq);
    }
}

void DepthBlock::divide(int n)
{
    if (n > 1)
    {
        _dp = (int) round((float) _dp / (float) n);
        _gq = (int) round((float) _gq / (float) n);
        _dpf = (int) round((float) _dpf / (float) n);
    }
}

int DepthBlock::size() const
{
    return (_end - _start + 1);
}

int DepthBlock::intersect_size(int rid, int start, int end) const
{
    if (_rid != rid)
    {
        return (0);
    }

    if (_end < start || _start > end)
    {
        return (0);
    }

    return (min(end, _end) - max(_start, start) + 1);
}

int DepthBlock::intersect_size(const DepthBlock &db) const
{
    return (intersect_size(db._rid, db._start, db._end));
}

DepthBlock DepthBlock::intersect(const DepthBlock &db)
{
    return (intersect(db._rid, db._start, db._end));
}

DepthBlock DepthBlock::intersect(int rid, int start, int end)
{
    if (!intersect_size(rid, start, end))
    {
        std::cerr << "DepthBlock " << _rid << ":" << _start + 1 << "-" << _end + 1 << " does not contain region: "
                  << rid << ":" << start + 1 << "-" << end + 1 << std::endl;
        ggutils::die("DepthBlock: bad coordinates.");
    }

    return {_rid, max(_start, start), min(end, _end), _dp, _dpf, _gq, _ploidy};
}

int DepthBlock::get_ploidy() {return _ploidy;}
