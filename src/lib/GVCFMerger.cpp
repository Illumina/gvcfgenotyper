#include "GVCFMerger.hh"
#include "utils.hh"
#include <htslib/hts.h>
#include <htslib/vcf.h>

extern "C" {
      size_t hts_realloc_or_die(unsigned long, unsigned long, unsigned long, unsigned long, int, void**, char const*);
}

//#define DEBUG

GVCFMerger::~GVCFMerger()
{
//    bcf_destroy(_output_record);
    hts_close(_output_file);
    bcf_hdr_destroy(_output_header);
    free(_format_gt);
    free(_format_ad);
    free(_format_dp);
    free(_format_gq);
    free(_format_gqx);
    free(_format_dpf);
    free(_format_ps);
    free(_format_adf);
    free(_format_adr);
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
                       const int is_file /*= 0*/)
{
    _num_variants=0;
    _num_gvcfs = input_files.size();
    _readers.reserve(_num_gvcfs);
    std::cerr << "Input GVCFs:" << std::endl;
    for (size_t i = 0; i < _num_gvcfs; i++)
    {
        std::cerr << input_files[i] << std::endl;
        _readers.emplace_back(input_files[i], reference_genome, buffer_size, region, is_file);
    }
    assert(_readers.size() == _num_gvcfs);

    _output_file = hts_open(!output_filename.empty() ? output_filename.c_str() : "-", ("w" + output_mode).c_str());

    if (!_output_file)
    {
        die("problem opening output file: " + output_filename);
    }

    size_t n_allele = 2;
    _format_pl = nullptr;
    _format_gt = (int32_t *) malloc(2 * _num_gvcfs * sizeof(int32_t));
    _format_ad = (int32_t *) malloc(n_allele * _num_gvcfs * sizeof(int32_t));

    _format_adf = (int32_t *) malloc(n_allele * _num_gvcfs * sizeof(int32_t));
    _format_adr = (int32_t *) malloc(n_allele * _num_gvcfs * sizeof(int32_t));

    _info_adf = (int32_t *) malloc(n_allele * sizeof(int32_t));
    _info_adr = (int32_t *) malloc(n_allele * sizeof(int32_t));
    _info_ac = (int32_t *) malloc(n_allele * sizeof(int32_t));

    _format_dp = (int32_t *) malloc(_num_gvcfs * sizeof(int32_t));
    _format_gq = (int32_t *) malloc(_num_gvcfs * sizeof(int32_t));
    _format_gqx = (int32_t *) malloc(_num_gvcfs * sizeof(int32_t));
    _format_ps = (int32_t *) malloc(_num_gvcfs * sizeof(int32_t));
    _format_dpf = (int32_t *) malloc(_num_gvcfs * sizeof(int32_t));
    _num_pl=_num_variants=0;
    build_header();
    _record_collapser.init(_output_header);
    _output_record = bcf_init1();
}

int GVCFMerger::get_next_variant()
{
    assert(_readers.size() == _num_gvcfs);
    assert(!all_readers_empty());
    bcf1_t *min_rec = nullptr;
    //int min_index = -1;
    for (size_t i = 0; i < _num_gvcfs; i++)
    {
        bcf1_t *rec = _readers[i].front();
        if (rec != nullptr)
        {
            if (min_rec == nullptr || bcf1_less_than(rec, min_rec))
            {
                min_rec = rec;
               // min_index = i;
            }
        }
    }
    assert(min_rec != nullptr);

    _record_collapser.setPosition(min_rec->rid, min_rec->pos);
    for (auto it = _readers.begin(); it != _readers.end(); it++)
    {
        auto variants = it->get_all_variants_in_interval(min_rec->rid, min_rec->pos);
        for (auto rec = variants.first; rec != variants.second; rec++)
        {
            if(get_variant_rank(*rec) == get_variant_rank(min_rec))
            {
                _record_collapser.allele(*rec);
            }
        }
    }

    return(1);
}

bool GVCFMerger::all_readers_empty()
{
    for (size_t i = 0; i < _num_gvcfs; i++)
    {
        if (!_readers[i].empty())
        {
            return (false);
        }
    }
    return (true);
}

void GVCFMerger::set_output_buffers_to_missing(int num_alleles)
{
    int ploidy = 2;
    _num_pl = get_number_of_likelihoods(ploidy,num_alleles)* _num_gvcfs;
    _format_pl = (int32_t *) realloc(_format_pl, _num_pl * sizeof(int32_t));
    std::fill(_format_pl, _format_pl + _num_pl, bcf_int32_missing);
    _format_ad = (int32_t *) realloc(_format_ad, num_alleles * _num_gvcfs * sizeof(int32_t));

    _format_adf = (int32_t *) realloc(_format_adf, num_alleles * _num_gvcfs * sizeof(int32_t));
    _format_adr = (int32_t *) realloc(_format_adr, num_alleles * _num_gvcfs * sizeof(int32_t));

    _info_adf = (int32_t *) realloc(_info_adf, num_alleles * sizeof(int32_t));
    _info_adr = (int32_t *) realloc(_info_adr, num_alleles * sizeof(int32_t));
    _info_ac = (int32_t *) realloc(_info_ac, num_alleles * sizeof(int32_t));

    std::fill(_format_gt, _format_gt + 2 * _num_gvcfs, bcf_gt_missing);
    std::fill(_format_ad, _format_ad + num_alleles * _num_gvcfs, 0);

    std::fill(_format_adf, _format_adf + num_alleles * _num_gvcfs, 0);
    std::fill(_format_adr, _format_adr + num_alleles * _num_gvcfs, 0);

    std::fill(_info_adf, _info_adf + num_alleles, 0);
    std::fill(_info_adr, _info_adr + num_alleles, 0);
    std::fill(_info_ac, _info_ac + num_alleles, 0);

    std::fill(_format_dp, _format_dp + _num_gvcfs, bcf_int32_missing);
    std::fill(_format_dpf, _format_dpf + _num_gvcfs, bcf_int32_missing);
    std::fill(_format_gq, _format_gq + _num_gvcfs, bcf_int32_missing);
    std::fill(_format_gqx, _format_gqx + _num_gvcfs, bcf_int32_missing);
    std::fill(_format_ps, _format_ps + _num_gvcfs, bcf_int32_missing);
}

bcf1_t *GVCFMerger::next()
{
    if (all_readers_empty())
    {
        return (nullptr);
    }
#ifdef DEBUG
    bcf1_t *prev_rec = nullptr;
    if(_num_variants>0)//DEBUG
    {
        prev_rec = bcf_dup(_output_record);//debug
        bcf_unpack(prev_rec,BCF_UN_ALL);
    }
#endif
    bcf_clear(_output_record);
    DepthBlock homref_block;//working structure to store homref info.
    get_next_variant(); //stores all the alleles at the next position.
    bcf_update_id(_output_header, _output_record, ".");
    _record_collapser.collapse(_output_record);
    _output_record->qual = 0;
#ifdef DEBUG
    std::cerr << "GVCFMerger::next()" << std::endl;
    print_variant(_output_header, _output_record);
    if(prev_rec!=nullptr && !bcf1_greater_than(_output_record,prev_rec))//DEBUG
    {
        print_variant(prev_rec);
        print_variant(_output_record);
        std::cerr<<bcf1_greater_than(_output_record,prev_rec)<<std::endl;
        std::cerr<<bcf1_equal(_output_record,prev_rec)<<std::endl;
        std::cerr<<bcf1_less_than(_output_record,prev_rec)<<std::endl;
        std::cerr<<get_variant_rank(_output_record)<< " "<<get_variant_rank(prev_rec) <<std::endl;
        throw std::runtime_error("incorrect variant order detected");
    }
#endif
    //fill in the format information for every sample.
    set_output_buffers_to_missing(_output_record->n_allele);

    // count the number of written PS tags
    unsigned nps_written(0);
    int32_t mean_mq = 0;
    int32_t num_mq = 0;
    for (size_t i = 0; i < _num_gvcfs; i++)
    {
        auto sample_variants = _readers[i].get_all_variants_up_to(_record_collapser.get_max());
        const bcf_hdr_t *sample_header = _readers[i].get_header();
        if (sample_variants.first!=sample_variants.second)
        {//this sample has variants at this position, we need to populate its FORMAT field
            int dst_genotype_count = 0; //this tracks how many destination (haploid) genotypes have been filled.
            for (auto it = sample_variants.first; it != sample_variants.second; it++)
            {
                bcf1_t *sample_record = *it;
                if(!bcf_float_is_missing(sample_record->qual))
                {
                    _output_record->qual += sample_record->qual;
                }
                Genotype g(_readers[i].get_header(), sample_record);
                int allele = _record_collapser.allele(sample_record);
                for (int genotype_index = 0; genotype_index < g._ploidy; genotype_index++)
                {
                    if ((sample_variants.second - sample_variants.first)== 1)//there is only one variant at this position in this sample. simple copy.
                    {
                       assert(dst_genotype_count<=2);
                        _format_gt[2 * i + dst_genotype_count] = bcf_gt_allele(g._gt[genotype_index]) == 0 ? bcf_gt_unphased(0) : bcf_gt_unphased(allele);
                        dst_genotype_count++;
                    }
                    else //there are multiple variants at this position. we need to do some careful genotype counting.
                    {
                        if(bcf_gt_allele(g._gt[genotype_index]) == 1)
                        {
                            if(dst_genotype_count>=2)
                            {
                                std::cerr << "WARNING: had to drop an allele in sample "+std::to_string(i)+" due to conflicting genotype calls" <<std::endl;
                                print_variant(sample_header,sample_record);
                            }
                            else
                            {
                                _format_gt[2 * i + dst_genotype_count] = bcf_gt_unphased(allele);
                                dst_genotype_count++;
                            }
                        }
                    }
                    _format_gq[i] = g.get_gq();
                    _format_gqx[i] = g.get_gqx();
                    if(g.is_dp_missing())
                    {
                        g.setDepthFromAD();
                    }
                    _format_dp[i] = g.get_dp();
                    _format_dpf[i] = g.get_dpf();
                    // AD
                    _format_ad[i*_output_record->n_allele] = g.get_ad(0);
                    _format_ad[i*_output_record->n_allele+allele] = g.get_ad(1);
                    // ADF
                    _format_adf[i*_output_record->n_allele] = g.get_adf(0);
                    _format_adf[i*_output_record->n_allele+allele] = g.get_adf(1);
                    // ADR
                    _format_adr[i*_output_record->n_allele] = g.get_adr(0);
                    _format_adr[i*_output_record->n_allele+allele] = g.get_adr(1);
                }
                int32_t sample_mq = 0;
                int nval=1;
                int32_t* ptr = &sample_mq;
                if (bcf_get_info_int32(sample_header,sample_record,"MQ",&ptr,&nval) > 0) {
                    mean_mq += sample_mq;
                    ++num_mq;
                }
            }
        }
        else    //this sample does not have the variant, reconstruct the format fields from homref blocks
        {
            _readers[i].get_depth(_output_record->rid, _output_record->pos, get_end_of_variant(_output_record), homref_block);
            _format_dp[i] = homref_block._dp;
            _format_dpf[i] = homref_block._dpf;
            _format_gq[i] = homref_block._gq;
            // GQX is missing for HOM REF
            _format_ad[i * _output_record->n_allele] = homref_block._dp;
            if (homref_block._dp > 0)
            {
                _format_gt[2 * i] = _format_gt[2 * i + 1] = bcf_gt_unphased(0);
            }
        }
        _readers[i].flush_buffer(_record_collapser.get_max());
        if((_format_gt[2*i]==bcf_gt_missing) != (_format_gt[2*i+1]==bcf_gt_missing))
        {
//            print_variant(_output_header,_output_record);
//            std::cerr << i << " " <<  bcf_gt_allele(_format_gt[2*i]) << "/" << bcf_gt_allele(_format_gt[2*i+1]) << std::endl;
         //   die("bad genotypes");
        }
    }
    _num_variants++;
#ifdef DEBUG
    std::cerr << "BUFFER SIZES:";
    for (size_t i = 0; i < _num_gvcfs; i++)
    {
        std::cerr << " (" << _readers[i].get_num_variants() << "," << _readers[i].get_num_depth() << ")";
    }
    std::cerr << std::endl;
#endif

    assert(bcf_update_genotypes(_output_header, _output_record, _format_gt, _num_gvcfs * 2)==0);
    assert(bcf_update_format_int32(_output_header, _output_record, "GQ", _format_gq, _num_gvcfs)==0);
    //int ret1 = bcf_update_format_int32(_output_header, _output_record, "GQX", _format_gqx, _num_gvcfs);
    assert(bcf_update_format_int32(_output_header, _output_record, "GQX", _format_gqx, _num_gvcfs)==0);
    if (nps_written > 0)
    {
//        bcf_update_format_int32(_output_header, _output_record, "PS", _format_ps, _num_gvcfs);
    }
    bcf_update_format_int32(_output_header, _output_record, "DP", _format_dp, _num_gvcfs);
    bcf_update_format_int32(_output_header, _output_record, "DPF", _format_dpf, _num_gvcfs);
    bcf_update_format_int32(_output_header, _output_record, "AD", _format_ad, _num_gvcfs * _output_record->n_allele);
    bcf_update_format_int32(_output_header, _output_record, "ADF", _format_adf, _num_gvcfs * _output_record->n_allele);
    bcf_update_format_int32(_output_header, _output_record, "ADR", _format_adr, _num_gvcfs * _output_record->n_allele);
    bcf_update_format_int32(_output_header, _output_record, "PL", _format_pl, _num_pl);

    // Write INFO/MQ
    if (num_mq>0) {
        mean_mq /= num_mq;
        bcf_update_info_int32(_output_header,_output_record,"MQ",&mean_mq,1);
    }

    // Calculate AC/AN using htslib standard functions
    int ret = bcf_calc_ac(_output_header,_output_record,_info_ac,BCF_UN_FMT);
    if (ret)
    {
        // sum over all allele counts to get AN
        int an = 0;
        for (int i=0; i<_output_record->n_allele; i++) {
            an += _info_ac[i];
        }
        bcf_update_info_int32(_output_header, _output_record, "AN", &an, 1);
        bcf_update_info_int32(_output_header, _output_record, "AC", _info_ac+1, _output_record->n_allele-1);
    } 

    // Calculate INFO/ADF + INFO/ADR
    for (size_t i=0;i<(_num_gvcfs*_output_record->n_allele);i+=_output_record->n_allele) {
        for (size_t j=0;j<_output_record->n_allele;++j) { 
            _info_adf[j] += _format_adf[i+j];
            _info_adr[j] += _format_adr[i+j];
        }
    }
    bcf_update_info_int32(_output_header,_output_record,"ADF",_info_adf,_output_record->n_allele);
    bcf_update_info_int32(_output_header,_output_record,"ADR",_info_adr,_output_record->n_allele);

    return (_output_record);
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
            print_variant(_output_header, _output_record);
            throw std::runtime_error("GVCFMerger::write_vcf variants out of order");
        }

        last_pos = _output_record->pos;
        last_rid = _output_record->rid;
        bcf_write1(_output_file, _output_header, _output_record);
        num_written++;
    }
    assert(all_readers_empty());
    std::cerr << "Wrote " << num_written << " variants" << std::endl;
}

void GVCFMerger::build_header()
{
    _output_header = bcf_hdr_init("w");
    bool force_samples = false;
    int repeat_count = 0;
    for (size_t i = 0; i < _num_gvcfs; i++)
    {
        const bcf_hdr_t *hr = _readers[i].get_header();
        for (int j = 0; j < bcf_hdr_nsamples(hr); j++)
        {
            string sample_name = hr->samples[j];
            if (bcf_hdr_id2int(_output_header, BCF_DT_SAMPLE, sample_name.c_str()) != -1)
            {
                if (force_samples)
                {
                    cerr << "Warning duplicate sample found.\t" << sample_name;
                    sample_name += ":R" + to_string(static_cast<long long>(repeat_count++));
                    cerr << " -> " << sample_name << endl;
                } else
                {
                    die("duplicate sample names. use --force-samples if you want to merge anyway");
                }
            }
            bcf_hdr_add_sample(_output_header, sample_name.c_str());
        }
    }

    bcf_hdr_append(_output_header, "##source=gvcfmerge-v0.0.0");
    bcf_hdr_append(_output_header, "##INFO=<ID=GN,Number=G,Type=Integer,Description=\"count of each genotype.\">");
    bcf_hdr_append(_output_header,
                   "##INFO=<ID=AD,Number=R,Type=Integer,Description=\"sum of allele depths for ALL individuals\">");
    bcf_hdr_append(_output_header,
                   "##INFO=<ID=ADF,Number=R,Type=Integer,Description=\"Sum of allelic depth on forward strand for ALL individuals\">");
    bcf_hdr_append(_output_header,
                   "##INFO=<ID=ADR,Number=R,Type=Integer,Description=\"Sum of allelic depth on reverse strand for ALL individuals\">");
    bcf_hdr_append(_output_header,
                   "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"sum of depth  across all samples\">");
    bcf_hdr_append(_output_header,
                   "##INFO=<ID=DPF,Number=1,Type=Integer,Description=\"sum of basecalls filtered from input prior to site genotyping\">");
    bcf_hdr_append(_output_header, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">");
    bcf_hdr_append(_output_header,
                   "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">");
    bcf_hdr_append(_output_header, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    bcf_hdr_append(_output_header,
                   "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Filtered basecall depth used for site genotyping\">");
    bcf_hdr_append(_output_header,
                   "##FORMAT=<ID=DPF,Number=1,Type=Integer,Description=\"Basecalls filtered from input prior to site genotyping\">");
    bcf_hdr_append(_output_header,
                   "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed.\">");
    bcf_hdr_append(_output_header, "##FORMAT=<ID=ADF,Number=.,Type=Integer,Description=\"Allelic depths on the forward strand\"");
    bcf_hdr_append(_output_header, "##FORMAT=<ID=ADR,Number=.,Type=Integer,Description=\"Allelic depths on the reverse strand\"");
    bcf_hdr_append(_output_header, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
    bcf_hdr_append(_output_header,
                   "##FORMAT=<ID=FT,Number=A,Type=Integer,Description=\"variant was PASS filter in original sample gvcf\">");
    bcf_hdr_append(_output_header,
                   "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification.\">");
    bcf_hdr_append(_output_header, "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set identifier\">");
    bcf_hdr_append(_output_header, "##FORMAT=<ID=GQX,Number=1,Type=Integer,Description=\"Empirically calibrated genotype quality score for variant sites, otherwise minimum of {Genotype quality assuming variant position,Genotype quality assuming non-variant position}\">");


    copy_contigs(_readers[0].get_header(), _output_header);
    bcf_hdr_write(_output_file, _output_header);
}
