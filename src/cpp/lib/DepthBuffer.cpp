#include "GVCFReader.hh"

//interpolates depth for a given interval a<=x<b
//returns 0 on success and -1 if the buffer didnt contain the interval
int DepthBuffer::Interpolate(const int rid, const int start, const int stop, DepthBlock &db)
{
    db.SetToMissing();
    auto dp_ptr = _buffer.begin();
    while (dp_ptr != _buffer.end() && (dp_ptr->end() < start || dp_ptr->rid() < rid))
    {
        dp_ptr++;
    }

    if (dp_ptr == _buffer.end())
    {
        return (-1);
    }
    if (dp_ptr->IntersectSize(rid,start,stop) <= 0) {
        //std::cerr << "DEPTH BUFFER INTERPOLATE ERROR" << std::endl;
        //dump();
        //assert(0);
        return(-1);
    }
    db = dp_ptr->Intersect(rid, start, stop);
    dp_ptr++;
    while (dp_ptr != _buffer.end() && dp_ptr->IntersectSize(rid, start, stop) > 0)
    {
        db.Add(dp_ptr->Intersect(rid, start, stop));
        dp_ptr++;
    }
    return (0);
}

void DepthBuffer::dump() 
{
	for (auto db: _buffer) {
		std::cerr << db.rid() << " " << db.start() << " " << db.end() << " " << db.dp() << "\n";
	}
}

void DepthBuffer::push_back(const DepthBlock& db)
{
    //sanity check on value being pushed
    if (!(_buffer.empty() || db.rid() != _buffer.back().rid() || db.start() == (1 + _buffer.back().end()) ||
          db.start() == _buffer.back().end()))
    {
        //if (!_buffer.empty())
        //{
        //    std::cerr << _buffer.back().rid() << ":" << _buffer.back().start() + 1 << "-" << _buffer.back().end() + 1 << "   ->   ";
        //    std::cerr << db.rid() << ":" << db.start() + 1 << "-" << db.end() + 1 << std::endl;         
        //}
        //ggutils::warn("non-contiguous homozygous reference blocks. Is this an Illumina GVCF?");
    }

    if (_buffer.empty() || db.rid() > _buffer.back().rid() || db.start() > _buffer.back().end())
    {
        _buffer.push_back(std::move(db));
    } 
}

int DepthBuffer::FlushBuffer(const int rid, const int pos)
{
    int num_flushed = 0;
    while (!_buffer.empty() && _buffer.front().rid() < rid)
    {
        _buffer.pop_front();
        num_flushed++;
    }
    while (!_buffer.empty() && _buffer.front().end() < pos && _buffer.front().rid() == rid)
    {
        _buffer.pop_front();
        num_flushed++;
    }
    return (num_flushed);
}

int DepthBuffer::FlushBuffer()
{
    int num_flushed=_buffer.size();
    _buffer.clear();
    return num_flushed;
}

DepthBlock *DepthBuffer::Back()
{
    if (_buffer.empty())
    {
        return (nullptr);
    }
    return (&_buffer.back());
}

size_t DepthBuffer::Size()
{
    return (_buffer.size());
}

bool DepthBuffer::Empty() {
    return _buffer.empty();
}
