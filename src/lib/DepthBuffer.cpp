#include "GVCFReader.hh"

//interpolates depth for a given interval a<=x<b
//returns 0 on success and -1 if the buffer didnt contain the interval
int DepthBuffer::interpolate(const int rid, const int start, const int stop, DepthBlock &db)
{
    db.set_missing();
    auto dp_ptr = _buffer.begin();
    while (dp_ptr != _buffer.end() && (dp_ptr->end() < start || dp_ptr->rid() < rid))
    {
        dp_ptr++;
    }

    if (dp_ptr == _buffer.end())
    {
        return (-1);
    }

    db = dp_ptr->intersect(rid, start, stop);
    dp_ptr++;
    while (dp_ptr != _buffer.end() && dp_ptr->intersect_size(rid, start, stop) > 0)
    {
        db.add(dp_ptr->intersect(rid, start, stop));
        dp_ptr++;
    }
    return (0);
}

void DepthBuffer::push_back(const DepthBlock& db)
{
    //sanity check on value being pushed
    if (!(_buffer.empty() || db.rid() != _buffer.back().rid() || db.start() == (1 + _buffer.back().end()) ||
          db.start() == _buffer.back().end()))
    {
        if (!_buffer.empty())
        {
            std::cerr << db.rid() << ":" << db.start() + 1 << "-" << db.end() + 1 << "   ->   " << _buffer.back().rid()
                      << ":" << _buffer.back().start() + 1 << "-" << _buffer.back().end() + 1 << std::endl;
        }
        die("DepthBuffer: bad homref block");
    }

    if (_buffer.empty() || db.rid() > _buffer.back().rid() || db.start() > _buffer.back().end())
    {
        _buffer.push_back(std::move(db));
    }
}

int DepthBuffer::flush_buffer(const int rid, const int pos)
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

int DepthBuffer::flush_buffer()
{
    int num_flushed=_buffer.size();
    _buffer.clear();
    return num_flushed;
}

DepthBlock *DepthBuffer::back()
{
    if (_buffer.empty())
    {
        return (nullptr);
    }
    return (&_buffer.back());
}

size_t DepthBuffer::size()
{
    return (_buffer.size());
}
