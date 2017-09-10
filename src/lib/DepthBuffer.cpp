#include "GVCFReader.hpp"
#include "GVCFReader.hpp"

//interpolates depth for a given interval a<=x<b
int DepthBuffer::interpolate(int rid,int start,int stop,DepthBlock & db)
{
    db.zero();
    auto dp_ptr = _buffer.begin();
    while(dp_ptr != _buffer.end() && dp_ptr->_end < start)
    {
	dp_ptr++;
    }
    if(dp_ptr == _buffer.end())
    {
	die("dp buffer over run");
    }
    int num_intervals = 0;
    while(dp_ptr != _buffer.end() && dp_ptr->intersect_size(rid,start,stop)>0)
    {
	db.add(*dp_ptr);
	dp_ptr++;
	num_intervals++;
    }
    db.divide(num_intervals);
    return(0);
}

int DepthBuffer::push_back(DepthBlock db)
{
//    std::cerr << db._rid << ":" <<db._start+1<<"-"<<db._end+1<<std::endl;
    assert(_buffer.empty() || db._rid!=_buffer.back()._rid || db._start==(1+_buffer.back()._end));
    _buffer.push_back(db);
}

int DepthBuffer::flush_buffer(int rid,int pos)    
{
    while(_buffer.size()>0 &&  _buffer.front()._end <= pos && _buffer.front()._rid <= rid)
    {
	_buffer.pop_front();
    }    
}

int DepthBuffer::flush_buffer()
{
    flush_buffer(_buffer.back()._rid,_buffer.back()._end);
}
