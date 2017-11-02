#include "VariantBuffer.hh"

VariantBuffer::VariantBuffer()
{
    _num_duplicated_records = 0;
}

VariantBuffer::~VariantBuffer()
{
    flush_buffer();
}

bool VariantBuffer::has_variant(bcf1_t *v)
{
    if (_buffer.empty())
    {
        return (false);
    }
    int i = _buffer.size() - 1;
    while (i >= 0 && _buffer[i]->pos >= v->pos)
    {
        if (bcf1_equal(v, _buffer[i]))
        {
            return (true);
        }
        i--;
    }
    return (false);
}

int VariantBuffer::push_back(bcf1_t *rec)
{
    bcf_unpack(rec, BCF_UN_ALL);
    if (has_variant(rec))
    {
        _num_duplicated_records++;
        bcf_destroy(rec);
        return (0);
    }

    _buffer.push_back(rec);
    //moves the new record back through the buffer until buffer is sorted ie. one iteration of insert-sort
    int i = _buffer.size() - 1;
    while (i > 0 && bcf1_less_than(_buffer[i], _buffer[i - 1]))
    {
        bcf1_t *tmp = _buffer[i - 1];
        _buffer[i - 1] = _buffer[i];
        _buffer[i] = tmp;
        i--;
    }
    return (1);
}

int VariantBuffer::flush_buffer(bcf1_t *record)
{
    assert(record!=nullptr);
    int num_flushed = 0;
    while (!_buffer.empty() && bcf1_leq(_buffer.front(), record))
    {
        bcf_destroy(_buffer.front());
        _buffer.pop_front();
        num_flushed++;
    }
    return (num_flushed);
}

int VariantBuffer::flush_buffer(int chrom, int pos)
{
    int num_flushed = 0;
    while (!_buffer.empty() && _buffer.front()->rid < chrom)
    {
        bcf_destroy(_buffer.front());
        _buffer.pop_front();
        num_flushed++;
    }
    while (!_buffer.empty() && _buffer.front()->pos < pos && _buffer.front()->rid == chrom)
    {
        bcf_destroy(_buffer.front());
        _buffer.pop_front();
        num_flushed++;
    }
    return (num_flushed);
}

int VariantBuffer::flush_buffer()
{
    if (!_buffer.empty())
    {
        int ret = flush_buffer(_buffer.back()->rid, _buffer.back()->pos);
        return (ret);
    }
    else
    {
        return (0);
    }
}

size_t VariantBuffer::size()
{
    return (_buffer.size());
}

bool VariantBuffer::empty()
{
    return (_buffer.empty());
}


bcf1_t *VariantBuffer::back()
{
    if (_buffer.empty())
    {
        return (nullptr);
    }
    else
    {
        bcf1_t *ret = _buffer.back();
        assert(ret != nullptr);
        bcf_unpack(ret, BCF_UN_ALL);
        return (ret);
    }
}

bcf1_t *VariantBuffer::front()
{
    if (_buffer.empty())
    {
        return (nullptr);
    }
    else
    {
        bcf1_t *ret = _buffer.front();
        assert(ret != nullptr);
        bcf_unpack(ret, BCF_UN_ALL);
        return (ret);
    }
}


bcf1_t *VariantBuffer::pop()
{
    if (_buffer.empty())
    {
        return (nullptr);
    }
    else
    {
        bcf1_t *ret = _buffer.front();
        _buffer.pop_front();
        return (ret);
    }
}

pair<std::deque<bcf1_t *>::iterator,std::deque<bcf1_t *>::iterator> VariantBuffer::get_all_variants_in_interval(int chrom,int stop)
{
    auto a = _buffer.begin();
    pair<std::deque<bcf1_t *>::iterator,std::deque<bcf1_t *>::iterator > ret(a, a);
    if(_buffer.empty() || (*a)->rid!=chrom || stop < (*a)->pos )
    {
        return(ret);
    }

    while(ret.second!=_buffer.end() && chrom==(*ret.second)->rid && stop>=(*ret.second)->pos)
    {
        ret.second++;
    }
    return(ret);
}

pair<std::deque<bcf1_t *>::iterator,std::deque<bcf1_t *>::iterator> VariantBuffer::get_all_variants_up_to(bcf1_t * record)
{
    auto a = _buffer.begin();
    pair<std::deque<bcf1_t *>::iterator,std::deque<bcf1_t *>::iterator> ret(a, a);
    if(_buffer.empty() || bcf1_greater_than(*a,record) )
    {
        return(ret);
    }
    while(ret.second!=_buffer.end() && bcf1_leq(*ret.second,record))
    {
        ret.second++;
    }
    return(ret);
}
