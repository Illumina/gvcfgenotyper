#include "Genotype.hh"

#include <stdexcept>


int Genotype::get_gq()
{
    return(*_gq);
}

int Genotype::get_gqx()
{
    return(*_gqx);
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

int Genotype::get_adf(int index)
{
    assert(index<_num_adf);
    return(_adf[index]);
}

int Genotype::get_adr(int index)
{
    assert(index<_num_adr);
    return(_adr[index]);
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
    _num_adf = _num_allele;
    _num_adr = _num_allele;
    _num_gq = 1;
    _num_gqx = 1;
    _num_dp = 1;
    _num_dpf = 1;
    _gt = ggutils::zeros(_num_gt);
    _pl = ggutils::zeros(_num_pl);
    _ad = ggutils::zeros(_num_ad);
    _adf = ggutils::zeros(_num_ad);
    _adr = ggutils::zeros(_num_ad);
    _gq = ggutils::zeros(_num_gq);
    _gqx = ggutils::zeros(_num_gqx);
    _dp = ggutils::zeros(_num_dp);
    _dpf = ggutils::zeros(_num_dpf);
    _gl.assign(_num_pl, 0.);
    _adf_found=false;
    _adr_found=false;
}

Genotype::Genotype(bcf_hdr_t const *header, bcf1_t *record)
{
    _num_allele = record->n_allele;
    bcf_unpack(record, BCF_UN_ALL);
    assert(_num_allele > 1);

    //this chunk of codes reads our canonical FORMAT fields (PL,GQ,DP,DPF,AD)
    _gt = (int *)malloc(2*sizeof(int));//force gt to be of length 2.
    _gt[1] = bcf_int32_vector_end;
    _ad = _adf = _adr = _gq = _gqx = _dpf = _dp = _pl = nullptr;
    _num_gt = 2;
    _ploidy = bcf_get_genotypes(header, record, &_gt, &_num_gt);
    assert(_ploidy >= 0 && _ploidy <= 2);
    _num_pl = _num_gl = _ploidy == 1 ? _num_allele : _num_allele * (1 + _num_allele) / 2;
    _num_ad = 0, _num_adf = 0, _num_adr = 0, _num_gq = 0, _num_gqx = 0, _num_dpf = 0,  _num_dp = 0;

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
        ggutils::print_variant((bcf_hdr_t *)header,record);
        std::cerr << "Got " << ret << " values instead of " << _num_gl << " ploidy="<<_ploidy<<" num_allele="<<_num_allele<<std::endl;
        throw std::runtime_error("incorrect number of values in  FORMAT/PL");
    }
    assert(bcf_get_format_int32(header, record, "AD", &_ad, &_num_ad) == _num_allele);
    if(bcf_get_format_int32(header, record, "ADF", &_adf, &_num_adf) == _num_allele) {
        _adf_found=true;
    } else {
        std::cerr << "WARNING: gvcf without ADF supplied." << "\n";
        _adf_found=false;
    }
    if(bcf_get_format_int32(header, record, "ADR", &_adr, &_num_adr) == _num_allele) {
        _adr_found=true;
    } else {
        std::cerr << "WARNING: gvcf without ADR supplied." << "\n";
        _adr_found=false;
    }
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
    bcf_get_format_int32(header, record, "GQX", &_gqx, &_num_gqx);
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
        _gl[i] =ggutils::unphred(_pl[i]);
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
    memcpy(ret._gqx, _gqx, sizeof(int32_t) * _num_gqx);

    //marginalises FORMAT/AD
    ret._ad[reference_allele] = _ad[reference_allele];
    ret._ad[primary_allele] = _ad[index];
    ret._ad[symbolic_allele] = 0;

    if (_adf_found) {
        ret._adf[reference_allele] = _adf[reference_allele];
        ret._adf[primary_allele] = _adf[index];
        ret._adf[symbolic_allele] = 0;
    }

    if (_adr_found) {
        ret._adr[reference_allele] = _adr[reference_allele];
        ret._adr[primary_allele] = _adr[index];
        ret._adr[symbolic_allele] = 0;
    }
    for (int j = 1; j < _num_allele; j++)
    {
        if (index != j)//j is an alternate allele that is not i
        {
            ret._ad[symbolic_allele] += _ad[j];
            if (_adf_found) {
                ret._adf[symbolic_allele] += _adf[j];
            }
            if (_adr_found) {
                ret._adr[symbolic_allele] += _adr[j];
            }
        }
    }
    //FIXME: dummy values for the time being because PLs are complicated.
    for (int j = 0; j < ret._num_pl; j++)
    {
        ret._pl[j] =ggutils::phred(ret._gl[j]);
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
    free(_adf);
    free(_adr);
    free(_gq);
    free(_gqx);
    free(_dpf);
    free(_pl);
    free(_dp);
}

int Genotype::update_bcf1_t(bcf_hdr_t *header, bcf1_t *record)
{
    float max_gl = *std::max_element(_gl.begin(), _gl.end());
    for (int i = 0; i < _num_pl; i++)
    {
        _pl[i] = _gl[i] > 0 ? ggutils::phred(_gl[i] / max_gl) : 255;
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

    assert(bcf_update_format_int32(header, record, "GQX", _gqx, _num_gqx)==0);

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
    bcf_update_format_int32(header, record, "ADR", _adr, _num_adr);
    bcf_update_format_int32(header, record, "ADF", _adf, _num_adf);
    bcf_update_format_int32(header, record, "PL", _pl, _num_pl);

    return (0);
}

Genotype Genotype::collapse_alleles_into_ref(vector<int> & indices)
{
    int num_new_allele = indices.size()+1;
    for(size_t i=0;i<indices.size();i++)
    {
        assert(indices[i]>0 && indices[i]<_num_allele);
    }

    Genotype ret(_ploidy, num_new_allele);
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
    memcpy(ret._gqx, _gqx, sizeof(int32_t) * _num_gqx);

    std::vector<int> allele_map;
    allele_map.push_back(0);
    int dst_allele=1;
    for(int i=1;i<_num_allele;i++)
    {
        //if allele i is not our list, we collapse its values into the reference allele
        if (std::find(indices.begin(), indices.end(), i) == indices.end())
            allele_map.push_back(0);
        else
            allele_map.push_back(dst_allele++);
    }

    for(int i=0;i<_num_allele;i++)
    {
        ret._ad[allele_map[i]]+=_ad[i];
        ret._adr[allele_map[i]]+=_adf[i];
        ret._adf[allele_map[i]]+=_adr[i];
        for(int j=i;j<_num_allele;j++)
            ret._gl[ggutils::get_gl_index(allele_map[i],allele_map[j])] += _gl[ggutils::get_gl_index(i,j)];
    }
    ret.PLfromGL();

    bool is_phased = _ploidy==1 ? bcf_gt_is_phased(_gt[0]) : bcf_gt_is_phased(_gt[1]) && bcf_gt_is_phased(_gt[0]);

    for(int i=0;i<_ploidy;i++)
    {
        if(is_phased)
            ret._gt[i] = bcf_gt_phased(allele_map[bcf_gt_allele(_gt[i])]);
        else
            ret._gt[i] = bcf_gt_unphased(allele_map[bcf_gt_allele(_gt[i])]);
    }

    return(ret);
}

void Genotype::PLfromGL()
{
    float max_gl = *std::max_element(_gl.begin(), _gl.end());
    for (int i = 0; i < _num_pl; i++)
    {
        _pl[i] = _gl[i] > 0 ? ggutils::phred(_gl[i] / max_gl) : 255;
//        _pl[i] = _pl[i] > 255 ? 255 : _pl[i]; //sets a minimum of 255 on PL
    }
}
