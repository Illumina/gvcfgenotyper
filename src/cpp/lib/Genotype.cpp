#include "Genotype.hh"

#include <stdexcept>
#include <htslib/vcf.h>

#define MAXPL 255

void Genotype::SetDp(int val)
{
    _dp = (int32_t *)realloc(_dp,sizeof(int32_t));
    *_dp = val;
}

void Genotype::SetDpf(int val)
{
    _dpf = (int32_t *)realloc(_dpf,sizeof(int32_t));
    *_dpf = val;
}

void Genotype::SetGq(int val)
{
    _gq = (int32_t *)realloc(_gq,sizeof(int32_t));
    *_gq = val;
}

void Genotype::SetGqx(int val)
{
    _gqx = (int32_t *)realloc(_gqx,sizeof(int32_t));
    *_gqx = val;
}

int Genotype::gq()
{
    assert(_num_gq==1);
    return(*_gq);
}

int Genotype::gqx()
{
    assert(_num_gqx==1 || _num_gqx==0);
    if(_num_gqx==1)
        return(*_gqx);
    else
        return(bcf_int32_missing);
}

int Genotype::dp()
{
    if(_num_dp>0)
        return(*_dp);
    else
        return(bcf_int32_missing);
}

int Genotype::ad(int index)
{
    assert(index<_num_ad && index>=0);
    return(_ad[index]);
}

int Genotype::adf(int index)
{
    if(_num_adf==0) return(bcf_int32_missing);
    assert(index<_num_adf && index>=0);
    return(_adf[index]);
}

int Genotype::adr(int index)
{
    if(_num_adr==0) return(bcf_int32_missing);
    assert(index<_num_adr && index>=0);
    return(_adr[index]);
}

int Genotype::dpf()
{
    if(_num_dpf>0)
        return(*_dpf);
    else
        return(bcf_int32_missing);
}

bool Genotype::IsDpMissing()
{
    return(_num_dp==0);
}

Genotype::Genotype(int ploidy, int num_allele)
{
    allocate(ploidy,num_allele);
}

void Genotype::CallGenotype()
{
    if(_ploidy==1)
    {
        auto it = std::max_element(_gl.begin(), _gl.end());
        _gt[0] = bcf_gt_unphased(it - _gl.begin());
    }
    else
    {
        int maxi=0,maxj=0;
        for(int i=0;i<_num_allele;i++)
        {
            for(int j=i;j<_num_allele;j++)
            {
                std::cerr<<_gl[ggutils::get_gl_index(i,j)]<<std::endl;//debug
                if(_gl[ggutils::get_gl_index(i,j)]>_gl[ggutils::get_gl_index(maxi,maxj)])
                {
                    maxi=i;
                    maxj=j;
                }
            }
        }
        _gt[0] = bcf_gt_unphased(maxi);
        _gt[1] = bcf_gt_unphased(maxj);
    }
}

void Genotype::allocate(int ploidy, int num_allele)
{
    _has_pl=true;
    _ploidy = ploidy;
    _num_allele = num_allele;
    assert(_ploidy >= 0 && _ploidy <= 2);
    _num_gt = _ploidy;
    _num_pl = _ploidy == 1 ? _num_allele : _num_allele * (1 + _num_allele) / 2;
    _num_ad = _num_allele;
    _num_adf = _num_allele;
    _num_adr = _num_allele;
    _num_gq = 1;
    _num_gqx = 1;
    _num_dp = 1;
    _num_dpf = 1;
    _gt = (int32_t *)realloc(_gt,_ploidy*sizeof(int32_t));
    std::fill(_gt,_gt+_ploidy,bcf_gt_missing);
    _pl = ggutils::assign_bcf_int32_missing(_pl,_num_pl);
    _ad = ggutils::assign_bcf_int32_missing(_ad,_num_ad);
    _adf = ggutils::assign_bcf_int32_missing(_adf,_num_ad);
    _adr = ggutils::assign_bcf_int32_missing(_adr,_num_ad);
    _gq = ggutils::assign_bcf_int32_missing(_gq,_num_gq);
    _gqx = ggutils::assign_bcf_int32_missing(_gqx,_num_gqx);
    _dp = ggutils::assign_bcf_int32_missing(_dp,_num_dp);
    _dpf = ggutils::assign_bcf_int32_missing(_dpf,_num_dpf);
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
    _gt = (int *) malloc(2 * sizeof(int));//force gt to be of length 2.
    _gt[1] = bcf_int32_vector_end;
    _ad = _adf = _adr = _gq = _gqx = _dpf = _dp = _pl = nullptr;
    _num_gt = 2;
    _ploidy = bcf_get_genotypes(header, record, &_gt, &_num_gt);
    assert(_ploidy >= 0 && _ploidy <= 2);
    _num_pl = _ploidy == 1 ? _num_allele : _num_allele * (1 + _num_allele) / 2;

    //These values with be updated by their respective bcf_get_* calls.
    _num_ad = 0, _num_adf = 0, _num_adr = 0, _num_gq = 0, _num_gqx = 0, _num_dpf = 0, _num_dp = 0;

    _qual = record->qual;
    _pl = (int32_t *) malloc(sizeof(int32_t) * _num_pl);
    int status;
    status = bcf_get_format_int32(header, record, "PL", &_pl, &_num_pl);
    if (status == 1 || status == -3 || status == -1)
    {
        std::fill(_pl, _pl + _num_pl, bcf_int32_missing);
        status = _num_pl;
        _has_pl = false;
    }
    else
    {
        _has_pl=true;
    }

    if (status != _num_pl)
    {
        ggutils::print_variant((bcf_hdr_t *) header, record);
        std::cerr << "Got " << status << " values instead of " << _num_pl << " ploidy=" << _ploidy << " num_allele=" << _num_allele << std::endl;
        ggutils::die("incorrect number of values in  FORMAT/PL");
    }
    assert(bcf_get_format_int32(header, record, "AD", &_ad, &_num_ad) == _num_allele);
    if (bcf_get_format_int32(header, record, "ADF", &_adf, &_num_adf) == _num_allele)
    {
        _adf_found = true;
    }
    else
    {
        _adf_found = false;
    }
    if (bcf_get_format_int32(header, record, "ADR", &_adr, &_num_adr) == _num_allele)
    {
        _adr_found = true;
    }
    else
    {
        _adr_found = false;
    }
    status = bcf_get_format_int32(header, record, "DP", &_dp, &_num_dp);
    if (status != 1)
    {
        if (status == -3)
            SetDepthFromAd();
        else
            ggutils::die("problem extracting FORMAT/DP");
    }
    ggutils::bcf1_get_one_info_int(header,record,"MQ",_mq);
    status = bcf_get_format_int32(header, record, "DPF", &_dpf, &_num_dpf);
    if (status != 1)
    {
        _dpf = (int32_t *)malloc(sizeof(int32_t));
        *_dpf = bcf_int32_missing;
        _num_dpf=1;
    }
    bcf_get_format_int32(header, record, "GQX", &_gqx, &_num_gqx);

    if (bcf_get_format_int32(header, record, "GQ", &_gq, &_num_gq) == -2)
    {
        _gq = ggutils::assign_bcf_int32_missing(_gq,1);
        float *tmp_gq = nullptr;
        if(bcf_get_format_float(header, record, "GQ", &tmp_gq, &_num_gq) != 1)
        {
            //FIXME: we should be logging this in a separate file
            //std::cerr<<"WARNING: missing GQ value at pos "<<record->pos+1<<std::endl;
        }
        else
        {
            _gq[0] = (int32_t)tmp_gq[0];
            free(tmp_gq);
        }
        _num_gq = 1;
    }
    SetGlFromPl();
}

void Genotype::SetGlFromPl()
{
    _gl.assign(_num_pl, 1.);
    float den = 0.;
    for (int i = 0; i < _num_pl; i++)
    {
        _gl[i] =ggutils::unphred(_pl[i]);
        den += _gl[i];
    }
    for (int i = 0; i < _num_pl; i++) _gl[i] /= den;
}


void Genotype::SetDepthFromAd()
{
    if(_num_dp==0)
    {
        _dp = (int *)malloc(sizeof(int));
        _num_dp=1;
    }
    _dp[0] = 0;
    for (int i = 0; i < _num_ad; i++) _dp[0] += _ad[i];
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

int Genotype::UpdateBcfRecord(bcf_hdr_t *header, bcf1_t *record)
{
//    float max_gl = *std::max_element(_gl.begin(), _gl.end());
//    for (int i = 0; i < _num_pl; i++)
//    {
//        _pl[i] = _gl[i] > 0 ? ggutils::phred(_gl[i] / max_gl) : MAXPL;
//    }

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
        bcf_update_format_int32(header, record, "DP", _dp, _num_dp);
        bcf_update_format_int32(header, record, "DPF", _dpf, _num_dpf);
    }

    bcf_update_format_int32(header, record, "AD", _ad, _num_ad);
    bcf_update_format_int32(header, record, "ADR", _adr, _num_adr);
    bcf_update_format_int32(header, record, "ADF", _adf, _num_adf);
    bcf_update_format_int32(header, record, "PL", _pl, _num_pl);

    return (0);
}

void Genotype::SetDepthToZero()
{
    std::fill(_dp,_dp+_num_dp,0);
    std::fill(_ad,_ad+_num_ad,0);
    std::fill(_adf,_adf+_num_adf,0);
    std::fill(_adr,_adr+_num_adr,0);
}

void Genotype::CollapseAllelesIntoRef(vector<int> &indices, Genotype &out)
{
    int num_new_allele = indices.size()+1;
    for(size_t i=0;i<indices.size();i++)
    {
        assert(indices[i]>0 && indices[i]<_num_allele);
    }

    out.allocate(_ploidy, num_new_allele);
    out.SetDepthToZero();

    if(_num_dp>0)
        memcpy(out._dp, _dp, sizeof(int32_t) * _num_dp);
    if(_num_dpf>0)
        memcpy(out._dpf, _dpf, sizeof(int32_t) * _num_dpf);

    memcpy(out._gq, _gq, sizeof(int32_t) * _num_gq);
    memcpy(out._gqx, _gqx, sizeof(int32_t) * _num_gqx);

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
        out._ad[allele_map[i]]+=_ad[i];
        if(_adf_found)
            out._adf[allele_map[i]]+=_adf[i];
        if(_adr_found)
            out._adr[allele_map[i]]+=_adr[i];

        if(_ploidy==1)
            out._gl[allele_map[i]]+=_gl[i];
        else if(_ploidy==2)
            for(int j=i;j<_num_allele;j++)
                out._gl[ggutils::get_gl_index(allele_map[i],allele_map[j])] += _gl[ggutils::get_gl_index(i,j)];
        else
            ggutils::die("Genotype::CollapseAllelesIntoRef invalid ploidy");
    }
    out.SetPlFromGl();

    bool is_phased = _ploidy==1 ? bcf_gt_is_phased(_gt[0]) : bcf_gt_is_phased(_gt[1]) && bcf_gt_is_phased(_gt[0]);

    for(int i=0;i<_ploidy;i++)
    {
        if(bcf_gt_is_missing(_gt[i]))
            out._gt[i] = bcf_gt_missing;
        else if(is_phased)
            out._gt[i] = bcf_gt_phased(allele_map[bcf_gt_allele(_gt[i])]);
        else
            out._gt[i] = bcf_gt_unphased(allele_map[bcf_gt_allele(_gt[i])]);
    }
}

void Genotype::SetPlFromGl()
{
    float max_gl = *std::max_element(_gl.begin(), _gl.end());
    for (int i = 0; i < _num_pl; i++)
    {
        _pl[i] = _gl[i] > 0 ? ggutils::phred(_gl[i] / max_gl) : MAXPL; //fixes -0.0
        _pl[i] = _pl[i] > MAXPL ? MAXPL : _pl[i]; //sets a ceiling of MAXPL on PL
    }
}

int Genotype::PropagateFormatFields(size_t sample_index, size_t ploidy, ggutils::vcf_data_t *format)
{
    assert(sample_index<format->num_sample);
    //update scalars
    format->gq[sample_index] = gq();
    format->gqx[sample_index] = gqx();

    if(IsDpMissing()) SetDepthFromAd();

    format->dp[sample_index] = dp();
    format->dpf[sample_index] = dpf();

    //update Number=R (eg FORMAT/AD)
    for(int i=0;i<_num_allele;i++)
    {
        format->ad[sample_index * _num_allele + i] = ad(i);
        format->adf[sample_index * _num_allele + i] = adf(i);
        format->adr[sample_index * _num_allele + i] = adr(i);
    }
    format->gt[sample_index*ploidy] = gt(0);

    //update GT
    if(_ploidy==1)  format->gt[sample_index*ploidy+1] = bcf_int32_vector_end;
    else  format->gt[sample_index*ploidy+1] = gt(1);

    //update PL
    size_t num_pl_per_sample = ggutils::get_number_of_likelihoods(ploidy,_num_allele);
    size_t num_pl_in_this_sample = ggutils::get_number_of_likelihoods(_ploidy,_num_allele);
    int *dst = format->pl+num_pl_per_sample*sample_index;
    std::fill(dst,dst+num_pl_per_sample,bcf_int32_vector_end);
    memcpy(dst,_pl,num_pl_in_this_sample*sizeof(int));
    return(1);
}

void Genotype::SetGtToHomRef()
{
    _gt[0] = bcf_gt_unphased(0);
    if(_ploidy==2) _gt[1]=bcf_gt_unphased(0);
}

void Genotype::MakeDiploid()
{
    assert(_ploidy==1);
    _ploidy=2;
    _gt = (int32_t *)realloc(_gt,_ploidy*sizeof(int32_t));
    _gt[1] = bcf_gt_unphased(bcf_gt_allele(_gt[0]));
    int _new_num_pl = _ploidy == 1 ? _num_allele : _num_allele * (1 + _num_allele) / 2;
    int32_t *_new_pl = (int32_t *)malloc(sizeof(int32_t)*_new_num_pl);
    std::fill(_new_pl,_new_pl+_new_num_pl,MAXPL);
    for(int i=0;i< num_allele();i++)
        _new_pl[ggutils::get_gl_index(i,i)] = _pl[i];
    SetGlFromPl();
    free(_pl);
    _pl = _new_pl;
}

Genotype::Genotype(bcf_hdr_t *sample_header,
                   pair<std::deque<bcf1_t *>::iterator,std::deque<bcf1_t *>::iterator> & sample_variants,
                   multiAllele & alleles_to_map)
{
    int ploidy=0;
    size_t num_sample_variants = (sample_variants.second - sample_variants.first);
    for (auto it = sample_variants.first; it != sample_variants.second; it++) ploidy = max(ggutils::get_ploidy(sample_header,*it),ploidy);
    assert(ploidy==1 || ploidy==2);
    allocate(ploidy, alleles_to_map.GetNumAlleles()+1);
    SetDepthToZero();
    SetGtToHomRef();
    assert(_num_allele>1);
    std::fill(_pl,_pl+ggutils::get_number_of_likelihoods(ploidy,_num_allele),MAXPL);
    int dst_genotype_count=0;
    _qual = 0;
    _mq = 0;
    bool recall_genotypes_from_gl = false; //if this is triggered we well reset GT as argmax(GL). This is a hack to get around allele collisions (which are usually in flakey genomic regions anyway).
    for (auto it = sample_variants.first; it != sample_variants.second; it++)
    {
        Genotype g(sample_header, *it);
        if(g._ploidy!=ploidy)
        {
            std::cerr<<"WARNING: conflicting ploidy for sample "<<sample_header->samples[0]<<" "<< bcf_hdr_id2name(sample_header, (*it)->rid);
            std::cerr<< ":" << (*it)->pos + 1 << std::endl;
            if(g.ploidy()==1)
                g.MakeDiploid();
            else
                ggutils::die("Unresolvable ploidy conflict.");
        }

        //FIXME: These values are overwriting on each iteration. Ideally they should be the same so it does not matter,
        // but there will be situations where this is not the case due variants shifting position. Hence with have to
        // add some hacks to handle conflicting values.
        _qual = max(_qual, g.qual()); //FIXME: QUAL should be estimated from PL
        _mq = g.mq()==bcf_int32_missing ? bcf_int32_missing : max(mq(), g.mq());
        *_gq = *_gq==bcf_int32_missing ? g.gq() : max(gq(), g.gq());
        *_gqx = *_gqx==bcf_int32_missing ? g.gqx() : max(gqx(), g.gqx());
        *_dpf = *_dpf==bcf_int32_missing ? g.dpf() : *_dpf;

        _ad[0]  = max(ad(0), g.ad(0));
        _adf[0] = max(adf(0), g.adf(0));
        _adr[0] = max(adr(0), g.adr(0));
        _gl[0] = g.ploidy()==1 ? g.gl(0) : g.gl(0,0);
        _pl[0] = g.ploidy()==1 ? g.pl(0) : g.pl(0,0);

        //these are values specific to the alt allele in question.
        int dst_allele_index = alleles_to_map.Allele(*it);
        assert(dst_allele_index < _num_allele && dst_allele_index > 0);
        _ad[dst_allele_index] = g.ad(1);
        if(g.HasAdf() && g.HasAdr())
        {
            _adf[dst_allele_index] = g.adf(1);
            _adr[dst_allele_index] = g.adr(1);
        }
        if(g.HasPl())
        {
            if (ploidy == 1)
            {
                _gl[dst_allele_index] = g.gl(1);
            }
            else
            {
                for (auto allele2 = sample_variants.first; allele2 != sample_variants.second; allele2++)
                {
                    int dst_allele_index2 = alleles_to_map.Allele(*allele2);
                    int src_allele_index2 = ggutils::find_allele(*it,*allele2,1);
                    _gl[ggutils::get_gl_index(dst_allele_index2, dst_allele_index)] = g.gl(src_allele_index2, 1);
                    _pl[ggutils::get_gl_index(dst_allele_index2, dst_allele_index)] = g.pl(src_allele_index2, 1);
                }
                _gl[ggutils::get_gl_index(0, dst_allele_index)] = g.gl(0, 1);
                _pl[ggutils::get_gl_index(0, dst_allele_index)] = g.pl(0, 1);
            }
        }

        for (int genotype_index = 0; genotype_index < _ploidy; genotype_index++)
        {

            if (num_sample_variants == 1)//only one allele so this is a straightforward copy
            {
                if (bcf_gt_allele(g.gt(genotype_index)) == 0)
                    _gt[dst_genotype_count] = bcf_gt_unphased(0);
                if (bcf_gt_allele(g.gt(genotype_index)) == 1)
                    _gt[dst_genotype_count] = bcf_gt_unphased(dst_allele_index);
                dst_genotype_count++;
            }
            else //there are multiple variants at this position. we need to do some careful genotype counting.
            {
                if (bcf_gt_allele(g.gt(genotype_index)) == 1)
                {
                    if(dst_genotype_count>=_ploidy)
                    {
                        std::cerr<<"WARNING: conflicting alleles/genotypes for sample "<<sample_header->samples[0]<<" ";
                        std::cerr<< bcf_hdr_id2name(sample_header, (*it)->rid) << ":" << (*it)->pos + 1 << std::endl;
                        recall_genotypes_from_gl=true;
                    }
                    else
                    {
                        _gt[dst_genotype_count] = bcf_gt_unphased(dst_allele_index);
                        dst_genotype_count++;
                    }
                }
            }
        }
    }
    SetDepthFromAd();
    if(recall_genotypes_from_gl)
    {
        SetPlFromGl();
        CallGenotype();
    }

    //This is a sanity check which should possibly reside within unit testing.
    if(_ploidy==2 && (_gt[0]==bcf_gt_missing) != (_gt[1]==bcf_gt_missing))
    {
        ggutils::print_variant(sample_header,*sample_variants.first);
        ggutils::die("bad genotypes at sample "+(string)sample_header->samples[0]);
    }
}

float Genotype::qual()
{
    return(_qual);
}

int Genotype::mq()
{
    return(_mq);
}

int Genotype::ploidy() {return _ploidy;};
int Genotype::num_allele() {return _num_allele;};

int Genotype::gt(int index)
{
    assert(index<_ploidy);
    return(_gt[index]);
};

int Genotype::pl(int g0, int g1)
{
    assert(_ploidy==2);
    assert(g0>=0 && g1>=0);
    assert(g0<_num_allele && g1<_num_allele);
    return(_pl[ggutils::get_gl_index(g0,g1)]);
}

int Genotype::pl(int g0)
{
    assert(_ploidy==1);
    assert(g0>=0);
    assert(g0<_num_allele);
    return(_pl[g0]);
}

float Genotype::gl(int g0, int g1)
{
    if(_num_pl<=0) return(bcf_float_missing);
    assert(_ploidy==2);
    assert(g0>=0 && g1>=0);
    assert(g0<_num_allele && g1<_num_allele);
    return(_gl[ggutils::get_gl_index(g0,g1)]);
}

float Genotype::gl(int g0)
{
    if(_num_pl<=0) return(bcf_float_missing);
    assert(_ploidy==1);
    assert(g0>=0);
    assert(g0<_num_allele);
    return(_gl[g0]);
}
