#include <htslib/vcf.h>
#include "GVCFReader.hh"
#include "StringUtil.hh"
//#define DEBUG

int GVCFReader::flush_buffer(bcf1_t *record)
{
    assert(record!=nullptr);
    _depth_buffer.flush_buffer(record->rid, record->pos - 1);
    int num_flushed = _variant_buffer.flush_buffer(record);
    fill_buffer();
    return (num_flushed);
}

int GVCFReader::flush_buffer(int chrom, int pos)
{
    _depth_buffer.flush_buffer(chrom, pos);
    int num_flushed = _variant_buffer.flush_buffer(chrom, pos);
    fill_buffer();
    return (num_flushed);
}

int GVCFReader::flush_buffer()
{
    _depth_buffer.flush_buffer();
    return (_variant_buffer.flush_buffer());
}

GVCFReader::GVCFReader(const std::string &input_gvcf, const std::string &reference_genome_fasta, const int buffer_size,
                       const string &region /*=""*/, const int is_file /*=0*/)
{
    _bcf_record = nullptr;
    _bcf_reader = bcf_sr_init();
    if (!region.empty())
    {
        if (bcf_sr_set_regions(_bcf_reader, region.c_str(), is_file) == -1)
        {
            ggutils::die("Cannot navigate to region " + region);
        }
    }
    if (!(bcf_sr_add_reader(_bcf_reader, input_gvcf.c_str())))
    {
        ggutils::die("problem opening input");
    }
    if (buffer_size < 2)
    {
        ggutils::die("GVCFReader needs buffer size of at least 2");
    }
    _buffer_size = buffer_size;
    _bcf_record = nullptr;

    //header setup
    _bcf_header = _bcf_reader->readers[0].header;

    _normaliser = new Normaliser(reference_genome_fasta, _bcf_header);
    fill_buffer();

    // flush variant buffer to get rid of variants overlapping 
    // the interval start
    if(region.find(":")!=std::string::npos)
    {
        string chr;
        int64_t start, end = 0;
        stringutil::parsePos(region, chr, start, end);
        if (!region.empty())
        {
            int rid = bcf_hdr_name2id(_bcf_header, chr.c_str());
            flush_buffer(rid, start);
        }
    }
}

GVCFReader::~GVCFReader()
{
    if (_bcf_reader->errnum)
    {
        error("Error: %s\n", bcf_sr_strerror(_bcf_reader->errnum));
    }
    bcf_sr_destroy(_bcf_reader);
    delete _normaliser;
}

size_t GVCFReader::fill_buffer()
{

    if (_variant_buffer.size() < 2)
    {
        read_lines(2);
    }
    while (_variant_buffer.size() > 1 && _variant_buffer.back()->rid == _variant_buffer.front()->rid &&
           (_variant_buffer.back()->pos - _variant_buffer.front()->pos) < _buffer_size)
    {
        int num_read = read_lines(1);
        if (num_read == 0)
        {
            break;
        }
    }

    return (_variant_buffer.size());
}

int GVCFReader::read_until(int rid, int pos)
{
    int num_read = 0;
//    bcf1_t *rec=_variant_buffer.back();
    DepthBlock *db = _depth_buffer.back();
    while (db == nullptr)
    {
        read_lines(1);
        db = _depth_buffer.back();
    }

    while (db != nullptr && bcf_sr_has_line(_bcf_reader, 0) && (db->_rid < rid || (db->_rid == rid && db->_end < pos)))
    {
        if (read_lines(1) < 1)
        {
            break;
        }
        else
        {
            num_read++;
        }
        db = _depth_buffer.back();
    }

    return (num_read);
}

int GVCFReader::read_lines(const unsigned num_lines)
{
    if (num_lines == 0)
    {
        return (0);
    }

    unsigned num_read = 0;

    while (num_read < num_lines && bcf_sr_next_line(_bcf_reader))
    {
        _bcf_record = bcf_sr_get_line(_bcf_reader, 0);
        if(_bcf_record->n_allele>1)
        {
            int32_t pass = bcf_has_filter(_bcf_header, _bcf_record, (char *) ".");
            bcf_update_format_int32(_bcf_header, _bcf_record, "FT", &pass, 1);
            bcf_update_filter(_bcf_header, _bcf_record, nullptr, 0);
            bcf_update_id(_bcf_header, _bcf_record, nullptr);
            vector<bcf1_t *> atomised_variants; 
            _normaliser->unarise(_bcf_record,atomised_variants);
            for (auto v = atomised_variants.begin();v!=atomised_variants.end();v++)
            {
                _variant_buffer.push_back(_bcf_header,*v);
            }
            num_read++;
        } 

        //buffer a depth block
        int num_format_values = 0;
        int32_t *value_pointer = nullptr;
        if (bcf_get_format_int32(_bcf_header, _bcf_record, "DP", &value_pointer, &num_format_values) == 1)
        {//if DP is present, this is either a snp or a homref block and we want to store it in depth buffer;
            int start = _bcf_record->pos;
            int32_t dp, dpf, gq, end;
            dp = *value_pointer;
            end = ggutils::get_end_of_gvcf_block(_bcf_header, _bcf_record);

            //if it is a variant use GQ else use GQX (this is a illumina GVCF quirk)
            if (_bcf_record->n_allele>1)
            {
                bcf_get_format_int32(_bcf_header, _bcf_record, "GQ", &value_pointer, &num_format_values);//FIXME: we need code to handle float/int here. this is a general problem that we should consider using templating to solve
            }
            else
            {
                bcf_get_format_int32(_bcf_header, _bcf_record, "GQX", &value_pointer, &num_format_values);
            }
            gq = *value_pointer;
            bcf_get_format_int32(_bcf_header, _bcf_record, "DPF", &value_pointer, &num_format_values);
            dpf = *value_pointer;
            _depth_buffer.push_back(DepthBlock(_bcf_record->rid, start, end, dp, dpf, gq));
            free(value_pointer);
        }

    }
    return (num_read);
}

bcf1_t *GVCFReader::front()
{
    fill_buffer();
    return (_variant_buffer.front());
}

bcf1_t *GVCFReader::pop()
{
    int num_read = fill_buffer();

    if (_variant_buffer.empty() && num_read == 0)
    {
        return (nullptr);
    }

    return (_variant_buffer.pop());
}

bool GVCFReader::empty()
{
    if (_variant_buffer.empty())
    {
        return (fill_buffer() == 0);
    }
    else
    {
        return (false);
    }
}

const  bcf_hdr_t *GVCFReader::get_header()
{
    return (_bcf_header);
}

//gets dp/dpf/gq (possibly interpolated) for a give interval a<=x<b
void GVCFReader::get_depth(int rid, int start, int stop, DepthBlock &db)
{
    read_until(rid, stop);
    if (_depth_buffer.interpolate(rid, start, stop, db) < 0)
    {
        ggutils::die("GVCFReader::get_depth problem with depth buffer");
    }
}

size_t GVCFReader::get_num_variants()
{
    return (_variant_buffer.size());
}

size_t GVCFReader::get_num_depth()
{
    return (_depth_buffer.size());
}

//gets all variants in interval start<=x<=stop
pair<std::deque<bcf1_t *>::iterator,std::deque<bcf1_t *>::iterator> GVCFReader::get_all_variants_in_interval(int chrom,int stop)
{
    return(_variant_buffer.get_all_variants_in_interval(chrom,stop));
}

pair<std::deque<bcf1_t *>::iterator,std::deque<bcf1_t *>::iterator>  GVCFReader::get_all_variants_up_to(bcf1_t *record)
{
    return(_variant_buffer.get_all_variants_up_to(record));
}

