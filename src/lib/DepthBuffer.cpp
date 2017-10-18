#include "GVCFReader.hpp"

//interpolates depth for a given interval a<=x<b
//returns 0 on success and -1 if the buffer didnt contain the interval
int DepthBuffer::interpolate(int rid, int start, int stop, DepthBlock &db)
{
    db.set_missing();
    auto dp_ptr = _buffer.begin();
    while (dp_ptr != _buffer.end() && (dp_ptr->_end < start || dp_ptr->_rid < rid))
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

void DepthBuffer::push_back(DepthBlock db)
{
    //sanity check on value being pushed
    if (!(_buffer.empty() || db._rid != _buffer.back()._rid || db._start == (1 + _buffer.back()._end) ||
          db._start == _buffer.back()._end))
    {
        if (!_buffer.empty())
        {
            std::cerr << db._rid << ":" << db._start + 1 << "-" << db._end + 1 << "   ->   " << _buffer.back()._rid
                      << ":" << _buffer.back()._start + 1 << "-" << _buffer.back()._end + 1 << std::endl;
        }
        die("DepthBuffer: bad homref block");
    }

    if (_buffer.empty() || db._rid > _buffer.back()._rid || db._start > _buffer.back()._end)
    {
        _buffer.push_back(db);
    }
}

int DepthBuffer::flush_buffer(int rid, int pos)
{
    int num_flushed = 0;
    while (_buffer.size() > 0 && _buffer.front()._rid < rid)
    {
        _buffer.pop_front();
        num_flushed++;
    }
    while (_buffer.size() > 0 && _buffer.front()._end < pos && _buffer.front()._rid == rid)
    {
        _buffer.pop_front();
        num_flushed++;
    }
    return (num_flushed);
}

int DepthBuffer::flush_buffer()
{
    return (flush_buffer(_buffer.back()._rid, _buffer.back()._end));
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
