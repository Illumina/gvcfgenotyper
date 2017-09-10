#include "GVCFReader.hpp"

VariantBuffer::VariantBuffer() 
{
    _num_duplicated_records=0;
}

VariantBuffer::~VariantBuffer() 
{
    flush_buffer();
}

bool VariantBuffer::has_variant(bcf1_t *v) 
{
    int i = _buffer.size()-1;
    while(i>=0 && _buffer[i]->pos >= v->pos)
    {
	if(v==_buffer[i])
	{
	    return(true);
	}
	i--;
    }
    return(false);
}

int VariantBuffer::push_back(bcf1_t *v) 
{
    bcf_unpack(v, BCF_UN_ALL);
    if(has_variant(v))
    {
	_num_duplicated_records++;
	return(0);
    }

    _buffer.push_back(bcf_dup(v));
    int i = _buffer.size()-1;
    while(i>0 && _buffer[i]->pos < _buffer[i-1]->pos) 
    {
	bcf1_t *tmp=_buffer[i-1];
	_buffer[i-1]=_buffer[i];
	_buffer[i]=tmp;
	i--;
    }
    return(1);
}

int VariantBuffer::flush_buffer(int chrom,int pos)
{
    int num_flushed=0;
    while(_buffer.size()>0 &&  _buffer.front()->pos <= pos && _buffer.front()->rid <= chrom)
    {
	bcf1_t *rec = _buffer.front();	    
	_buffer.pop_front();
	num_flushed++;
    }    
    return(num_flushed);
}  

int VariantBuffer::flush_buffer()
{
    if(_buffer.size()>0) 
    {
	int ret = flush_buffer(_buffer.back()->rid,_buffer.back()->pos);
	return(ret);
    }
    else
    {
	return(0);
    }
}

size_t VariantBuffer::size()
{
    return(_buffer.size());
}

bool VariantBuffer::empty()
{
    return(_buffer.empty());
}

bcf1_t *VariantBuffer::front()
{
    if(_buffer.empty())
    {
	return(NULL);
    }
    else
    {	
	bcf1_t *ret = _buffer.front();
	bcf_unpack(ret, BCF_UN_ALL);
	return(ret);
    }
}


bcf1_t *VariantBuffer::pop()
{
    if(_buffer.empty())
    {
	return(NULL);
    }
    else
    {
	bcf1_t *ret = _buffer.front();
	bcf_unpack(ret, BCF_UN_ALL);
	_buffer.pop_front();
	return(ret);
    }
}
