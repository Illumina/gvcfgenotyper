#include "GVCFReader.hpp"

Genotype::Genotype(int ploidy, int num_allele)
{
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
    _num_gl = _ploidy == 1 ? _ploidy : _num_allele * (1 + _num_allele) / 2;
    _num_ad = 0, _num_gq = 0, _num_dpf = 0, _num_pl = 0, _num_dp = 0;

    int ret;
    ret = bcf_get_format_int32(header, record, "PL", &_pl, &_num_pl);
    assert(ret == _num_gl);
    assert(bcf_get_format_int32(header, record, "AD", &_ad, &_num_ad) == _num_allele);
    if (bcf_get_format_int32(header, record, "DP", &_dp, &_num_dp) == -3)
    {
        _dp = (int *) malloc(sizeof(int));
        _dp[0] = std::accumulate(_ad, _ad + _num_ad, 0);
    }
    if (bcf_get_format_int32(header, record, "DPF", &_dpf, &_num_dpf) == -3)
    {
        _dpf = (int *) malloc(sizeof(int));
        _dp[0] = 0;
    }

    if (bcf_get_format_int32(header, record, "GQ", &_gq, &_num_gq) ==
        -2)//catches cases where GQ is set as a float (should be int according to vcf4.3)
    {
        float *tmp_gq = nullptr;
        _num_gq = 0;
        assert(bcf_get_format_float(header, record, "GQ", &tmp_gq, &_num_gq) == 1);
        _gq = (int32_t  *)malloc(sizeof(int32_t));
        _gq[0] = (int32_t) tmp_gq[0];
        free(tmp_gq);
    }
    _gl.assign(_num_pl, 1.);
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

    memcpy(ret._dp, _dp, sizeof(int32_t) * _num_dp);
    memcpy(ret._dpf, _dpf, sizeof(int32_t) * _num_dpf);
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

    bcf_update_genotypes(header, record, _gt, _num_gt);
    bcf_update_format_int32(header, record, "GQ", _gq, _num_gq);
    bcf_update_format_int32(header, record, "DP", _dp, _num_dp);
    bcf_update_format_int32(header, record, "DPF", _dpf, _num_dpf);
    bcf_update_format_int32(header, record, "AD", _ad, _num_ad);
    bcf_update_format_int32(header, record, "PL", _pl, _num_pl);

    return (0);
}
