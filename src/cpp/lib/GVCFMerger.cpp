#include "GVCFMerger.hh"
#include "ggutils.hh"
#include <htslib/hts.h>
#include <htslib/vcf.h>

#include <unordered_map>

extern "C" {
      size_t hts_realloc_or_die(unsigned long, unsigned long, unsigned long, unsigned long, int, void**, char const*);
}

//#define DEBUG

GVCFMerger::~GVCFMerger()
{
    delete _normaliser;
    hts_close(_output_file);
    bcf_hdr_destroy(_output_header);
    delete _format;
    free(_info_adf);
    free(_info_adr);
    free(_info_ac);
    bcf_destroy(_output_record);
}

GVCFMerger::GVCFMerger(const vector<string> &input_files,
                       const string &output_filename,
                       const string &output_mode,
                       const string &reference_genome,
                       int buffer_size,
                       const string &region /*= ""*/,
                       const int is_file /*= 0*/,
                       bool ignore_non_matching_ref,
                       bool force_samples)
{    
    _force_samples = force_samples;
    _has_pl = true;
    _has_strand_ad=true;
    _num_variants=0;
    _normaliser = new Normaliser(reference_genome,ignore_non_matching_ref);
    _num_gvcfs = input_files.size();
    _readers.reserve(_num_gvcfs);

    // retrieve logger from factory
    _lg = spdlog::get("gg_logger");
    assert(_lg!=nullptr);
    _lg->info("Input GVCFs:");
    for (size_t i = 0; i < _num_gvcfs; i++)
    {
        _lg->info("Opened {} {}/{}",input_files[i],(i+1),_num_gvcfs);
        _readers.emplace_back(input_files[i], _normaliser, buffer_size, region, is_file);
        _has_pl &= _readers.back().HasPl();
        _has_strand_ad &= _readers.back().HasStrandAd();
    }
    assert(_readers.size() == _num_gvcfs);

    _output_file = hts_open(!output_filename.empty() ? output_filename.c_str() : "-", ("w" + output_mode).c_str());

    if (!_output_file)
    {
        ggutils::die("problem opening output file: " + output_filename);
    }

    size_t n_allele = 2;
    _format = new ggutils::vcf_data_t(2,2,_num_gvcfs);

    _info_adf = (int32_t *) malloc(n_allele * sizeof(int32_t));
    _info_adr = (int32_t *) malloc(n_allele * sizeof(int32_t));
    _info_ac = (int32_t *) malloc(n_allele * sizeof(int32_t));

    BuildHeader();
    _record_collapser.Init(_output_header);
    _output_record = bcf_init1();

    _num_ps_written = 0;
    _mean_mq = 0;
    _num_mq = 0;
    _max_alleles = INT32_MAX;
}

int GVCFMerger::GetNextVariant()
{
    assert(_readers.size() == _num_gvcfs);
    assert(!AreAllReadersEmpty());
    bcf1_t *min_rec = nullptr;
    //int min_index = -1;
    for (size_t i = 0; i < _num_gvcfs; i++)
    {
        bcf1_t *rec = _readers[i].Front();
        if (rec != nullptr)
        {
            if (min_rec == nullptr || ggutils::bcf1_less_than(rec, min_rec))
            {
                min_rec = rec;
               // min_index = i;
            }
        }
    }
    assert(min_rec != nullptr);

    _record_collapser.SetPosition(min_rec->rid, min_rec->pos);
    for(auto it = _readers.begin(); it != _readers.end(); it++)
    {
        auto variants = it->GetAllVariantsInInterval(min_rec->rid, min_rec->pos);
        for (auto rec = variants.first; rec != variants.second; rec++)
        {
            if(ggutils::get_variant_rank(*rec) == ggutils::get_variant_rank(min_rec))
            {
                _record_collapser.Allele(*rec);
            }
        }
    }
    assert(_record_collapser.GetNumAlleles()>0);
    return(1);
}

bool GVCFMerger::AreAllReadersEmpty()
{
    for (size_t i = 0; i < _num_gvcfs; i++)
    {
        if (!_readers[i].IsEmpty())
        {
            return (false);
        }
    }
    return (true);
}

void GVCFMerger::SetOutputBuffersToMissing(int num_alleles)
{
    _format->resize(num_alleles);
    _format->set_missing();

    _info_adf = (int32_t *) realloc(_info_adf, num_alleles * sizeof(int32_t));
    _info_adr = (int32_t *) realloc(_info_adr, num_alleles * sizeof(int32_t));
    _info_ac = (int32_t *) realloc(_info_ac, num_alleles * sizeof(int32_t));

    std::fill(_info_adf, _info_adf + num_alleles, 0);
    std::fill(_info_adr, _info_adr + num_alleles, 0);
    std::fill(_info_ac, _info_ac + num_alleles, 0);

}

void GVCFMerger::GenotypeHomrefVariant(int sample_index, DepthBlock &homref_block)
{
    int num_pl_per_sample = ggutils::get_number_of_likelihoods(2,_output_record->n_allele);
    int num_pl_in_this_sample = ggutils::get_number_of_likelihoods(homref_block.ploidy(),_output_record->n_allele);

    int *pl_ptr = _format->pl + sample_index*num_pl_per_sample;
    _format->dp[sample_index] = homref_block.dp();
    _format->dpf[sample_index] = homref_block.dpf();
    _format->gq[sample_index] = homref_block.gq();
    _format->ad[sample_index * _output_record->n_allele] = homref_block.dp();
    for(int i=1;i<_output_record->n_allele;i++)
        _format->ad[sample_index * _output_record->n_allele+i] = 0;
    if (homref_block.dp() > 0)
    {
        if(homref_block.ploidy()==2)
        {
            _format->gt[2 * sample_index] = _format->gt[2 * sample_index + 1] = bcf_gt_unphased(0);
        }
        else
        {
            _format->gt[2 * sample_index] = bcf_gt_unphased(0);
            _format->gt[2 * sample_index + 1] = bcf_int32_vector_end;
        }
    }
    //FIXME: dummy PL value for homref sites
    std::fill(pl_ptr,pl_ptr+num_pl_per_sample,bcf_int32_vector_end);
    std::fill(pl_ptr,pl_ptr+num_pl_in_this_sample,255);
    pl_ptr[0] = 0;
}

void GVCFMerger::GenotypeAltVariant(int sample_index,bcf1_t *sample_variants)
{
    int default_ploidy=2;
    Genotype g(_readers[sample_index].GetHeader(), sample_variants,_record_collapser);
    g.PropagateFormatFields(sample_index, default_ploidy, _format);
    if(g.mq() != bcf_int32_missing)
    {
        _mean_mq += g.mq();
        _num_mq++;
    }
    if(!bcf_float_is_missing(g.qual()))
        _output_record->qual += g.qual();
}

void GVCFMerger::GenotypeSample(int sample_index)
{
    DepthBlock homref_block;//working structure to store homref info.
    auto hdr = _readers[sample_index].GetHeader();
    auto records = _readers[sample_index].GetAllVariantsUpTo(_record_collapser.GetMax());

    bcf1_t *sample_record = CollapseRecords(hdr,records);
    //this sample has variants at this position, we need to populate its FORMAT field
    if (sample_record!=nullptr)
    {
        GenotypeAltVariant(sample_index, sample_record);
        bcf_destroy(sample_record);
    }
    else    //this sample does not have the variant, reconstruct the format fields from homref blocks
    {
        _readers[sample_index].GetDepth(_output_record->rid, _output_record->pos,
                                        ggutils::get_end_of_variant(_output_record), homref_block);
        GenotypeHomrefVariant(sample_index, homref_block);
    }
}

bcf1_t *GVCFMerger::next()
{
    if (AreAllReadersEmpty()) return (nullptr);

    bcf_clear(_output_record);
    GetNextVariant(); //stores all the alleles at the next position.
    bcf_update_id(_output_header, _output_record, ".");
    _record_collapser.Collapse(_output_record);
    _output_record->qual = 0;
#ifdef DEBUG
    ggutils::print_variant(_output_header,_output_record);
#endif

    if(_output_record->n_allele <= _max_alleles)
    {
        //fill in the format information for every sample.
        SetOutputBuffersToMissing(_output_record->n_allele);

        // count the number of written PS tags
        _num_ps_written = 0;
        _mean_mq = 0;
        _num_mq = 0;

        for (size_t i = 0; i < _num_gvcfs; i++)
        {
            GenotypeSample(i);
            _readers[i].FlushBuffer(_record_collapser.GetMax());
        }
        _num_variants++;
        UpdateFormatAndInfo();
        return(_output_record);
    }
    else
    {
        _lg->warn("Too many alleles at {}:{} dropping this position.",
                  bcf_hdr_id2name(_output_header,_output_record->rid),_output_record->pos+1);
        for (size_t i = 0; i < _num_gvcfs; i++)
            _readers[i].FlushBuffer(_record_collapser.GetMax());
        return(next());
    }
}

void GVCFMerger::setMedianInfoValues()
{
    std::vector<float> median_gq(_output_record->n_allele);
    std::vector<float> median_gqx(_output_record->n_allele);
    for(int allele=0;allele<_output_record->n_allele;allele++)
    {
        std::vector<int> index_of_alt_genotypes;
        for(size_t i=0;i<_num_gvcfs;i++)
        {
            bool is_alt;
            if(_format->ploidy==1)
                is_alt = !bcf_gt_is_missing(_format->gt[i]) && bcf_gt_allele(_format->gt[i])==allele;
            else
                is_alt = !bcf_gt_is_missing(_format->gt[2*i+1]) && !bcf_gt_is_missing(_format->gt[2*i]) && (bcf_gt_allele(_format->gt[2*i])==allele||bcf_gt_allele(_format->gt[2*i+1])==allele);
            if(is_alt)
                index_of_alt_genotypes.push_back(i);
        }
        std::vector<int> values_at_alt_genotypes;
        for(auto it=index_of_alt_genotypes.begin();it!=index_of_alt_genotypes.end();it++) 
            if(_format->gq[*it]!=bcf_int32_missing)
                values_at_alt_genotypes.push_back(_format->gq[*it]);
        bcf_float_set_missing(median_gq[allele]);
        if(!values_at_alt_genotypes.empty())
            median_gq[allele] =  ggutils::inplace_median(values_at_alt_genotypes);
        values_at_alt_genotypes.clear();
        for(auto it=index_of_alt_genotypes.begin();it!=index_of_alt_genotypes.end();it++) 
            if(_format->gqx[*it]!=bcf_int32_missing)
                values_at_alt_genotypes.push_back(_format->gqx[*it]);

        bcf_float_set_missing(median_gqx[allele]);
        if(!values_at_alt_genotypes.empty())
            median_gqx[allele] =  ggutils::inplace_median(values_at_alt_genotypes);
    }
    assert(bcf_update_info_float(_output_header,_output_record,"GQX_MEDIAN",median_gqx.data()+1,_output_record->n_allele-1)==0);
    assert(bcf_update_info_float(_output_header,_output_record,"GQ_MEDIAN",median_gq.data()+1,_output_record->n_allele-1)==0);
}

void GVCFMerger::UpdateFormatAndInfo()
{
    assert(bcf_update_genotypes(_output_header, _output_record,_format->gt, _num_gvcfs * 2)==0);
    assert(bcf_update_format_string(_output_header, _output_record, "FT",(const char **)_format->ft, _num_gvcfs)==0);    
    assert(bcf_update_format_int32(_output_header, _output_record, "GQ",_format->gq, _num_gvcfs)==0);
    assert(bcf_update_format_int32(_output_header, _output_record, "GQX",_format->gqx, _num_gvcfs)==0);
    assert(bcf_update_format_int32(_output_header, _output_record, "DP",_format->dp, _num_gvcfs)==0);
    assert(bcf_update_format_int32(_output_header, _output_record, "DPF",_format->dpf, _num_gvcfs)==0);
    assert(bcf_update_format_int32(_output_header, _output_record, "AD",_format->ad, _num_gvcfs * _output_record->n_allele)==0);
    if(_has_strand_ad)
    {
        assert(bcf_update_format_int32(_output_header, _output_record, "ADF",_format->adf, _num_gvcfs * _output_record->n_allele)==0);
        assert(bcf_update_format_int32(_output_header, _output_record, "ADR",_format->adr, _num_gvcfs * _output_record->n_allele)==0);
    }
    if(_has_pl) assert(bcf_update_format_int32(_output_header, _output_record, "PL",_format->pl, _format->num_pl)==0);

    // Write INFO/MQ
    if (_num_mq>0)
    {
        _mean_mq /= _num_mq;
        assert(bcf_update_info_int32(_output_header,_output_record,"MQ",&_mean_mq,1)==0);
    }

    // Calculate AC/AN using htslib standard functions
    int ret = bcf_calc_ac(_output_header,_output_record,_info_ac,BCF_UN_FMT);
    if (ret)
    {
        // sum over all allele counts to get AN
        int an = 0;
        for (int i=0; i<_output_record->n_allele; i++)
            an += _info_ac[i];

        bcf_update_info_int32(_output_header, _output_record, "AN", &an, 1);
        bcf_update_info_int32(_output_header, _output_record, "AC", _info_ac+1, _output_record->n_allele-1);
    }

    // Calculate INFO/ADF + INFO/ADR
    if(_has_strand_ad)
    {
        for (size_t i=0;i<(_num_gvcfs*_output_record->n_allele);i+=_output_record->n_allele)
        {
            for (size_t j=0;j<_output_record->n_allele;++j)
            {
                assert( (_format->adr[i+j]==bcf_int32_missing) == (_format->adf[i+j]==bcf_int32_missing) );
                if(_format->adf[i+j]!=bcf_int32_missing)
                {
                    _info_adf[j] += _format->adf[i+j];
                    _info_adr[j] += _format->adr[i+j];
                }
            }
        }
        bcf_update_info_int32(_output_header,_output_record,"ADF",_info_adf,_output_record->n_allele);
        bcf_update_info_int32(_output_header,_output_record,"ADR",_info_adr,_output_record->n_allele);
        ggutils::fisher_sb_test(_info_adr,_info_adf,_output_record->n_allele,_sb_pvalue);
        bcf_update_info_float(_output_header,_output_record,"FS",_sb_pvalue.data(),_output_record->n_allele-1);
    }
    setMedianInfoValues();
}

void GVCFMerger::write_vcf()
{
    int last_rid = -1;
    int last_pos = 0;
    int num_written = 0;
    while (next() != nullptr)
    {
        if (!(_output_record->pos >= last_pos || _output_record->rid > last_rid))
        {
            std::cerr << bcf_hdr_int2id(_output_header, BCF_DT_CTG, _output_record->rid) << ":"
                      << _output_record->pos + 1 << std::endl;
            std::cerr << bcf_hdr_int2id(_output_header, BCF_DT_CTG, last_rid) << ":" << last_pos + 1 << std::endl;
            ggutils::print_variant(_output_header, _output_record);
            throw std::runtime_error("GVCFMerger::write_vcf variants out of order");
        }

        last_pos = _output_record->pos;
        last_rid = _output_record->rid;
        bcf_write1(_output_file, _output_header, _output_record);
        num_written++;
    }
    assert(AreAllReadersEmpty());
    _lg->info("Wrote {} variants",num_written);
}

void GVCFMerger::BuildHeader()
{
    _output_header = bcf_hdr_init("w");
    std::unordered_map<std::string,long long> repeat_count;
    for (size_t i = 0; i < _num_gvcfs; i++)
    {
        const bcf_hdr_t *hr = _readers[i].GetHeader();
        for (int j = 0; j < bcf_hdr_nsamples(hr); j++)
        {
            string sample_name = hr->samples[j];
            if (bcf_hdr_id2int(_output_header, BCF_DT_SAMPLE, sample_name.c_str()) != -1)
            {
                if (_force_samples)
                {
                    _lg->warn("Warning duplicate sample found.\t {} -> {}",
                              sample_name,
                              (sample_name + ":R" + to_string(repeat_count[sample_name]))
                              ); 
                    // Need to modify sample name here and not inside the log message
                    sample_name += ":R" + to_string(repeat_count[sample_name]++);
                } else
                {
                    ggutils::die("duplicate sample names. use --force-samples if you want to merge anyway");
                }
            }
            bcf_hdr_add_sample(_output_header, sample_name.c_str());
        }
    }
    
    bcf_hdr_append(_output_header, "##INFO=<ID=FS,Number=A,Type=Float,Description=\"Fisher exact test for per allele strand bias.\">");
    bcf_hdr_append(_output_header, "##INFO=<ID=GN,Number=G,Type=Integer,Description=\"count of each genotype.\">");
    bcf_hdr_append(_output_header,
                   "##INFO=<ID=AD,Number=R,Type=Integer,Description=\"sum of allele depths for ALL individuals\">");
    bcf_hdr_append(_output_header,
                   "##INFO=<ID=ADF,Number=R,Type=Integer,Description=\"Sum of allelic depth on forward strand for ALL individuals\">");
    bcf_hdr_append(_output_header,
                   "##INFO=<ID=ADR,Number=R,Type=Integer,Description=\"Sum of allelic depth on reverse strand for ALL individuals\">");
    bcf_hdr_append(_output_header,
                   "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"sum of depth  across all samples\">");
    bcf_hdr_append(_output_header,"##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"RMS of mapping quality\">");
    bcf_hdr_append(_output_header,
                   "##INFO=<ID=DPF,Number=1,Type=Integer,Description=\"sum of basecalls filtered from input prior to site genotyping\">");
    bcf_hdr_append(_output_header, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">");
    bcf_hdr_append(_output_header,
                   "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">");
    bcf_hdr_append(_output_header,
                   "##INFO=<ID=GQX_MEDIAN,Number=A,Type=Float,Description=\"The median GQX value for samples carrying this alternate allele\">");
    bcf_hdr_append(_output_header,
                   "##INFO=<ID=GQ_MEDIAN,Number=A,Type=Float,Description=\"The median GQ value for samples carrying this alternate allele\">");
                   
    bcf_hdr_append(_output_header, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    bcf_hdr_append(_output_header,
                   "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Filtered basecall depth used for site genotyping\">");
    bcf_hdr_append(_output_header,
                   "##FORMAT=<ID=DPF,Number=1,Type=Integer,Description=\"Basecalls filtered from input prior to site genotyping\">");
    bcf_hdr_append(_output_header,
                   "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed.\">");
    bcf_hdr_append(_output_header, "##FORMAT=<ID=ADF,Number=R,Type=Integer,Description=\"Allelic depths on the forward strand\"");
    bcf_hdr_append(_output_header, "##FORMAT=<ID=ADR,Number=R,Type=Integer,Description=\"Allelic depths on the reverse strand\"");
    bcf_hdr_append(_output_header, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
    bcf_hdr_append(_output_header, "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Sample filter, 'PASS' indicates that all single sample filters passed for this sample\">");

    bcf_hdr_append(_output_header,
                   "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in "
                           "the VCF specification.\">");
    bcf_hdr_append(_output_header, "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set identifier\">");
    bcf_hdr_append(_output_header, "##FORMAT=<ID=GQX,Number=1,Type=Integer,Description=\"Empirically calibrated genotype quality score for "
            "variant sites, otherwise minimum of {Genotype quality assuming variant position,Genotype quality assuming non-variant position}\">");
    bcf_hdr_append(_output_header, ("##gvcfgenotyper_version="+(string)GG_VERSION).c_str());
    ggutils::copy_contigs(_readers[0].GetHeader(), _output_header);
    
    bcf_hdr_write(_output_file, _output_header);
}
