#include "GVCFMerger.hpp"
#include "utils.hpp"
#include <stdexcept>
// #define DEBUG

GVCFMerger::~GVCFMerger()
{
//    bcf_destroy(_output_record);
    hts_close(_output_file);
    bcf_hdr_destroy(_output_header);
    free(_format_gt);
    free(_format_ad);
    free(_format_dp);
    free(_format_gq);
    free(_format_dpf);
    free(_format_ps);
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
    _num_gvcfs = input_files.size();
    _readers.reserve(_num_gvcfs);
    std::cerr << "Input GVCFs:" << std::endl;
    for (int i = 0; i < _num_gvcfs; i++)
    {
        std::cerr << input_files[i] << std::endl;
        _readers.emplace_back(input_files[i], reference_genome, buffer_size, region, is_file);
    }
    assert(_readers.size() == _num_gvcfs);

    _output_file = hts_open(output_filename != "" ? output_filename.c_str() : "-", ("w" + output_mode).c_str());

    if (!_output_file)
    {
        die("problem opening output file: " + output_filename);
    }

    _format_gt = (int32_t *) malloc(2 * _num_gvcfs * sizeof(int32_t));
    _format_ad = (int32_t *) malloc(2 * _num_gvcfs * sizeof(int32_t));
    _format_dp = (int32_t *) malloc(_num_gvcfs * sizeof(int32_t));
    _format_gq = (int32_t *) malloc(_num_gvcfs * sizeof(int32_t));
    _format_ps = (int32_t *) malloc(_num_gvcfs * sizeof(int32_t));
    _format_dpf = (int32_t *) malloc(_num_gvcfs * sizeof(int32_t));

    build_header();
    _output_record = bcf_init1();
}

bcf1_t *GVCFMerger::get_next_variant()
{
    assert(_readers.size() == _num_gvcfs);
    if (all_readers_empty())
    {
        return (nullptr);
    }
    bcf1_t *min_rec = nullptr;
    int min_index = -1;
    for (int i = 0; i < _num_gvcfs; i++)
    {
        bcf1_t *rec = _readers[i].front();
        if (rec != nullptr)
        {
            if (min_rec == nullptr || bcf1_less_than(rec, min_rec))
            {
                min_rec = rec;
                min_index = i;
            }
        }
    }

    assert(min_rec != nullptr);
    return (min_rec);
}

bool GVCFMerger::all_readers_empty()
{
    for (int i = 0; i < _num_gvcfs; i++)
    {
        if (!_readers[i].empty())
        {
            return (false);
        }
    }
    return (true);
}

void GVCFMerger::set_output_buffers_to_missing()
{
    for (int i = 0; i < _num_gvcfs; i++)
    {
        _format_gt[i * 2] = bcf_gt_missing;
        _format_gt[i * 2 + 1] = bcf_gt_missing;
        _format_ad[i * 2] = bcf_int32_missing;
        _format_ad[i * 2 + 1] = bcf_int32_missing;
        _format_dp[i] = bcf_int32_missing;
        _format_gq[i] = bcf_int32_missing;
        _format_dpf[i] = bcf_int32_missing;
        _format_ps[i] = bcf_int32_missing;
    }
}

bcf1_t *GVCFMerger::next()
{
    if (all_readers_empty())
    {
        return (nullptr);
    }
    bcf_clear(_output_record);
    int n_allele = 2;
    DepthBlock homref_block;//working structure to store homref info.
    //copy allele information from new variant
    bcf1_t *next_variant = get_next_variant();
    assert(bcf1_not_equal(next_variant, _output_record));
    assert(next_variant != nullptr);
    bcf_update_id(_output_header, _output_record, ".");
    _output_record->rid = next_variant->rid;
    _output_record->pos = next_variant->pos;
    bcf_update_alleles(_output_header, _output_record, (const char **) next_variant->d.allele, next_variant->n_allele);
    _output_record->qual = 0;
#ifdef DEBUG
    print_variant(_output_header,_output_record);
#endif
    //fill in the format information for every sample.
    set_output_buffers_to_missing();

    // count the number of written PS tags
    unsigned nps_written(0);
    for (int i = 0; i < _num_gvcfs; i++)
    {
        bcf1_t *sample_record = _readers[i].front();
        const bcf_hdr_t *sample_header = _readers[i].get_header();

        if (sample_record != nullptr && bcf1_equal(sample_record, _output_record))
        {//this sample has an explicit copy of the variant. just copy the format fields into output rrecorc
            _output_record->qual += sample_record->qual;
            int nval = 2;
            int32_t *ptr = _format_gt + 2 * i;
            int num_gt_in_sample = bcf_get_genotypes(sample_header, sample_record, &ptr, &nval);
            if (num_gt_in_sample != 2 && num_gt_in_sample != 1)
            {
                std::cerr << _output_record->rid << ":" << _output_record->pos + 1 << ":" << _output_record->d.allele[0]
                          << ":" << _output_record->d.allele[1] << std::endl;
                throw std::runtime_error("GVCFMerger: invalid genotypes");
            }
            if (num_gt_in_sample == 1)//pad the GT with vector_end
            {
                ptr[1] = bcf_int32_vector_end;
            }
            ptr = _format_ad + 2 * i;

            assert(bcf_get_format_int32(sample_header, sample_record, "AD", &ptr, &nval) == 2);
            nval = 1;
            ptr = _format_dp + i;
            if (bcf_get_format_int32(sample_header, sample_record, "DP", &ptr, &nval) < 0)
            {//FORMAT/DP not present (indels)P. we take DP = SUM(AD)
                _format_dp[i] = _format_ad[2 * i] + _format_ad[2 * i + 1];
            }
            ptr = _format_dpf + i;
            bcf_get_format_int32(sample_header, sample_record, "DPF", &ptr, &nval);

            //this is a bit dangerous: we specify GQ as Float when it should be Integer according to vcf spec
            //so we have to do some fiddly casting
            //solution is to make this general such that it handles integer or float
            float *gq_ptr = nullptr;
            nval = 0;
            if (bcf_get_format_float(sample_header, sample_record, "GQ", &gq_ptr, &nval) != 1)
            {
                std::cerr << "WARNING: missing FORMAT/GQ at " << bcf_hdr_id2name(_output_header, _output_record->rid) \
 << ":" << _output_record->pos + 1 << ":" << _output_record->d.allele[0] << ":" << _output_record->d.allele[1] \
 << std::endl;
            }
            else
            {
                _format_gq[i] = (int32_t) (*gq_ptr);
                free(gq_ptr);
            }
            ptr = _format_ps + i;
            int res = bcf_get_format_int32(sample_header, sample_record, "PS", &ptr, &nval);
            nps_written += (res > 0 ? res : 0);
        }
        else    //this sample does not have the variant, reconstruct the format fields from homref blocks
        {
            _readers[i].get_depth(_output_record->rid, _output_record->pos, get_end_of_variant(_output_record),
                                  homref_block);
            _format_dp[i] = homref_block._dp;
            _format_dpf[i] = homref_block._dpf;
            _format_gq[i] = homref_block._gq;
            _format_ad[i * 2] = homref_block._dp;
            _format_ad[i * 2 + 1] = 0;
            if (homref_block._dp > 0)
            {
                _format_gt[2 * i] = _format_gt[2 * i + 1] = bcf_gt_unphased(0);
            }
        }
        _readers[i].flush_buffer(_output_record);
    }
#ifdef DEBUG
    std::cerr << "BUFFER SIZES:";
    for(int i=0;i<_num_gvcfs;i++)
    {
    std::cerr << " ("<<_readers[i].get_num_variants()<<","<<_readers[i].get_num_depth()<<")";
    }
    std::cerr<<std::endl;
#endif

    bcf_update_genotypes(_output_header, _output_record, _format_gt, _num_gvcfs * 2);
    bcf_update_format_int32(_output_header, _output_record, "GQ", _format_gq, _num_gvcfs);
    if (nps_written > 0)
    {
        bcf_update_format_int32(_output_header, _output_record, "PS", _format_ps, _num_gvcfs);
    }
    bcf_update_format_int32(_output_header, _output_record, "DP", _format_dp, _num_gvcfs);
    bcf_update_format_int32(_output_header, _output_record, "DPF", _format_dpf, _num_gvcfs);
    bcf_update_format_int32(_output_header, _output_record, "AD", _format_ad, _num_gvcfs * n_allele);
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
    for (int i = 0; i < _num_gvcfs; i++)
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
                }
                else
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
    bcf_hdr_append(_output_header, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
    bcf_hdr_append(_output_header,
                   "##FORMAT=<ID=FT,Number=A,Type=Integer,Description=\"variant was PASS filter in original sample gvcf\">");
    bcf_hdr_append(_output_header,
                   "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification.\">");
    bcf_hdr_append(_output_header, "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set identifier\">");

    copy_contigs(_readers[0].get_header(), _output_header);
    bcf_hdr_write(_output_file, _output_header);
}
