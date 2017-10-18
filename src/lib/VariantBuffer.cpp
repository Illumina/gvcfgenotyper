#include "GVCFReader.hpp"

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

int VariantBuffer::flush_buffer(const bcf1_t *record)
{
    int num_flushed = 0;
    while (_buffer.size() > 0 && bcf1_leq(_buffer.front(), record))
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
    while (_buffer.size() > 0 && _buffer.front()->rid < chrom)
    {
        bcf_destroy(_buffer.front());
        _buffer.pop_front();
        num_flushed++;
    }
    while (_buffer.size() > 0 && _buffer.front()->pos < pos && _buffer.front()->rid == chrom)
    {
        bcf_destroy(_buffer.front());
        _buffer.pop_front();
        num_flushed++;
    }
    return (num_flushed);
}

int VariantBuffer::flush_buffer()
{
    if (_buffer.size() > 0)
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
