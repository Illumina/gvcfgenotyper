#include "Genotype.hh"

#include <stdexcept>


int Genotype::get_gq()
{
    return(*_gq);
}

int Genotype::get_dp()
{
    if(_num_dp>0)
    {
        return(*_dp);
    }
    else
    {
        return(bcf_int32_missing);
    }
}

int Genotype::get_ad(int index)
{
    assert(index<_num_ad);
    return(_ad[index]);
}

int Genotype::get_dpf()
{
    if(_num_dpf>0)
    {
        return(*_dpf);
    }
    else
    {
        return(bcf_int32_missing);
    }
}

void Genotype::set_dp_missing()
{
    _num_dp=0;
    _num_dpf=0;
}

bool Genotype::is_dp_missing()
{
    return(_num_dp==0);
}

Genotype::Genotype(int ploidy, int num_allele)
{
    _has_pl=true;
    _ploidy = ploidy;
    _num_allele = num_allele;
    assert(_ploidy >= 0 && _ploidy <= 2);
    _num_gt = _ploidy;
    _num_pl = _ploidy == 1 ? _ploidy : _num_allele * (1 + _num_allele) / 2;
    _num_ad = _num_allele;
    _num_gq = 1;
    _num_dp = 1;
    _num_dpf = 1;
    _gt = zeros(_num_gt);
    _pl = zeros(_num_pl);
    _ad = zeros(_num_ad);
    _gq = zeros(_num_gq);
    _dp = zeros(_num_dp);
    _dpf = zeros(_num_dpf);
    _gl.assign(_num_pl, 0.);
}

Genotype::Genotype(bcf_hdr_t const *header, bcf1_t *record)
{
    _num_allele = record->n_allele;
    bcf_unpack(record, BCF_UN_ALL);
    assert(_num_allele > 1);

    //this chunk of codes reads our canonical FORMAT fields (PL,GQ,DP,DPF,AD)
    _gt = (int *)malloc(2*sizeof(int));//force gt to be of length 2.
    _gt[1] = bcf_int32_vector_end;
    _ad = _gq = _dpf = _dp = _pl = nullptr;
    _num_gt = 2;
    _ploidy = bcf_get_genotypes(header, record, &_gt, &_num_gt);
    assert(_ploidy >= 0 && _ploidy <= 2);
    _num_pl = _num_gl = _ploidy == 1 ? _num_allele : _num_allele * (1 + _num_allele) / 2;
    _num_ad = 0, _num_gq = 0, _num_dpf = 0,  _num_dp = 0;

    _pl = (int32_t *)malloc(sizeof(int32_t)*_num_gl);
    int ret;
    ret = bcf_get_format_int32(header, record, "PL", &_pl, &_num_pl);
    if(ret==1 || ret==-3)
    {
        //std::cerr << "WARNING: missing FORMAT/PL at " << record->pos+1 <<std::endl;
        //std::cerr<<_pl[0]<<std::endl;
        //assert(_pl[0]==bcf_int32_missing);//FIXME: I don't understand why some times this is *not* bcf_int32_missing. Probably need to look at the htslib code to understand.
        std::fill(_pl,_pl+_num_gl,bcf_int32_missing);
        _num_pl=_num_gl;
        ret=_num_gl;
        _has_pl=false;
    }

    if(ret != _num_gl)
    {
        print_variant((bcf_hdr_t *)header,record);
        std::cerr << "Got " << ret << " values instead of " << _num_gl << " ploidy="<<_ploidy<<" num_allele="<<_num_allele<<std::endl;
        throw std::runtime_error("incorrect number of values in  FORMAT/PL");
    }
    assert(bcf_get_format_int32(header, record, "AD", &_ad, &_num_ad) == _num_allele);
    ret = bcf_get_format_int32(header, record, "DP", &_dp, &_num_dp);
    if(ret!=1)
    {
        if(ret==-3)
        {
            setDepthFromAD();
        }
        else
        {
            throw std::runtime_error("problem extracting FORMAT/DP");
        }
    }
    bcf_get_format_int32(header, record, "DPF", &_dpf, &_num_dpf);
    //FIXME: GQ is should be an integer but is sometimes set as float. we need to catch and handle this.
    if (bcf_get_format_int32(header, record, "GQ", &_gq, &_num_gq) == -2)
    {
        _gq = (int *)malloc(sizeof(int));
        float tmp_gq = -1;
        float *tmp_ptr = &tmp_gq;
        _num_gq = 1;
        if(bcf_get_format_float(header, record, "GQ", &tmp_ptr, &_num_gq) != 1)
        {
            std::cerr<<"WARNING: missing GQ value at pos "<<record->pos+1<<std::endl;
            _gq[0] = bcf_int32_missing;
        }
        else
        {
            _gq[0] = (int32_t)tmp_gq;
        }
        _num_gq = 1;
    }
    _gl.assign(_num_gl, 1.);
    for (int i = 0; i < _num_pl; i++)
    {
        _gl[i] = unphred(_pl[i]);
    }
}

Genotype Genotype::marginalise(int index)
{
    const int reference_allele = 0;
    const int primary_allele = 1;
    const int symbolic_allele = 2;
    const int num_new_allele = 3;
    assert(index > 0 && index < _num_allele);
    Genotype ret(_ploidy, num_new_allele);

    for (int j = 0; j < _ploidy; j++)
    {
        if (bcf_gt_allele(_gt[j]) == reference_allele)
        {
            ret._gt[j] = bcf_gt_unphased(reference_allele);
        }
        else if (bcf_gt_allele(_gt[j]) == index)
        {
            ret._gt[j] = bcf_gt_unphased(primary_allele);
        }
        else
        {
            ret._gt[j] = bcf_gt_unphased(symbolic_allele);
        }
    }
    if(_num_dp>0)
    {
        assert(_num_dpf>0);
        memcpy(ret._dp, _dp, sizeof(int32_t) * _num_dp);
        memcpy(ret._dpf, _dpf, sizeof(int32_t) * _num_dpf);
    }
    else
    {
        ret.set_dp_missing();
    }

    memcpy(ret._gq, _gq, sizeof(int32_t) * _num_gq);

    //marginalises FORMAT/AD
    ret._ad[reference_allele] = _ad[reference_allele];
    ret._ad[primary_allele] = _ad[index];
    ret._ad[symbolic_allele] = 0;
    for (int j = 1; j < _num_allele; j++)
    {
        if (index != j)//j is an alternate allele that is not i
        {
            ret._ad[symbolic_allele] += _ad[j];
        }
    }
    //FIXME: dummy values for the time being because PLs are complicated.
    for (int j = 0; j < ret._num_pl; j++)
    {
        ret._pl[j] = phred(ret._gl[j]);
    }

    return (ret);
}

void Genotype::setDepthFromAD()
{
    if(_num_dp==0)
    {
        _dp = (int *)malloc(sizeof(int));
        _dpf = (int *)malloc(sizeof(int));
        _num_dp=_num_dpf=1;
    }
    _dpf[0] = bcf_int32_missing;
    _dp[0] = 0;
    for (int i = 0; i < _num_ad; i++)
    {
        _dp[0] += _ad[i];
    }
}

Genotype::~Genotype()
{
    free(_gt);
    free(_ad);
    free(_gq);
    free(_dpf);
    free(_pl);
    free(_dp);
}

int Genotype::update_bcf1_t(bcf_hdr_t *header, bcf1_t *record)
{
    float max_gl = *std::max_element(_gl.begin(), _gl.end());
    for (int i = 0; i < _num_pl; i++)
    {
        _pl[i] = _gl[i] > 0 ? phred(_gl[i] / max_gl) : 255;
    }

    assert(bcf_update_genotypes(header, record, _gt, _num_gt)==0);

    //FIXME: hacking to resolve the float/int GQ issue
    if(bcf_hdr_id2type(header,BCF_HL_FMT,bcf_hdr_id2int(header,BCF_DT_ID,"GQ"))==BCF_HT_REAL)
    {
        float tmpgq =(float) _gq[0];
        int tmpnval = 1;
        assert(bcf_update_format_float(header, record, "GQ", &tmpgq, tmpnval)==0);
    }
    else
    {
        assert(bcf_update_format_int32(header, record, "GQ", _gq, _num_gq)==0);
    }

    if(_num_dp>0)
    {
        assert(_num_dpf>0);
        bcf_update_format_int32(header, record, "DP", _dp, _num_dp);
        bcf_update_format_int32(header, record, "DPF", _dpf, _num_dpf);
    }
    else
    {
        assert(_num_dpf==0);
        bcf_update_format_int32(header, record, "DP", NULL,0);
        bcf_update_format_int32(header, record, "DPF", NULL,0);
    }
    bcf_update_format_int32(header, record, "AD", _ad, _num_ad);
    bcf_update_format_int32(header, record, "PL", _pl, _num_pl);

    return (0);
}
