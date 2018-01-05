#include "VariantBuffer.hh"

VariantBuffer::VariantBuffer()
{
    _num_duplicated_records = 0;
}

VariantBuffer::~VariantBuffer()
{
    FlushBuffer();
}

bool VariantBuffer::HasVariant(const bcf_hdr_t *hdr, bcf1_t *v)
{
    if (_buffer.empty())
    {
        return (false);
    }
    int i = _buffer.size() - 1;
    while (i >= 0 && _buffer[i]->pos >= v->pos)
    {
        if (ggutils::bcf1_equal(v, _buffer[i]))
        {
            if (ggutils::is_hom_ref(hdr,_buffer[i]) && !ggutils::is_hom_ref(hdr,v)) {
                // if duplicate record in buffer is hom ref and the new record 
                // is not, swap them
                iter_swap(_buffer[i],v);
            }
            return (true);
        }
        i--;
    }
    return (false);
}

int VariantBuffer::PushBack(const bcf_hdr_t *hdr, bcf1_t *rec)
{
    bcf_unpack(rec, BCF_UN_ALL);
    if (HasVariant(hdr, rec))
    {
        _num_duplicated_records++;
        bcf_destroy(rec);
        return (0);
    }

    _buffer.push_back(rec);
    //moves the new record back through the buffer until buffer is sorted ie. one iteration of insert-sort
    int i = _buffer.size() - 1;
    while (i > 0 && ggutils::bcf1_less_than(_buffer[i], _buffer[i - 1]))
    {
        swap(_buffer[i - 1],_buffer[i]);
        i--;
    }
    return (1);
}

int VariantBuffer::FlushBuffer(bcf1_t *record)
{
    assert(record!=nullptr);
    int num_flushed = 0;
    while (!_buffer.empty() &&  ggutils::bcf1_leq(_buffer.front(), record))
    {
        bcf_destroy(_buffer.front());
        _buffer.pop_front();
        num_flushed++;
    }
    return (num_flushed);
}

int VariantBuffer::FlushBuffer(int chrom, int pos)
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

int VariantBuffer::FlushBuffer()
{
    int num_flushed = 0;
    while (!_buffer.empty())
    {
        bcf_destroy(_buffer.front());
        _buffer.pop_front();
        num_flushed++;
    }
    return num_flushed;
}

size_t VariantBuffer::Size()
{
    return (_buffer.size());
}

bool VariantBuffer::IsEmpty()
{
    return (_buffer.empty());
}


bcf1_t *VariantBuffer::Back()
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

bcf1_t *VariantBuffer::Front()
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


bcf1_t *VariantBuffer::Pop()
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

pair<std::deque<bcf1_t *>::iterator,std::deque<bcf1_t *>::iterator> VariantBuffer::GetAllVariantsInInterval(int chrom,
                                                                                                            int stop)
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

pair<std::deque<bcf1_t *>::iterator,std::deque<bcf1_t *>::iterator> VariantBuffer::GetAllVariantsUpTo(bcf1_t *record)
{
    auto a = _buffer.begin();
    pair<std::deque<bcf1_t *>::iterator,std::deque<bcf1_t *>::iterator> ret(a, a);
    if(_buffer.empty() ||  ggutils::bcf1_greater_than(*a,record) )
    {
        return(ret);
    }
    while(ret.second!=_buffer.end() &&  ggutils::bcf1_leq(*ret.second,record))
    {
        ret.second++;
    }
    return(ret);
}
