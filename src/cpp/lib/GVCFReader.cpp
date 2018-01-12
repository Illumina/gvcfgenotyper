#include <htslib/vcf.h>
#include "GVCFReader.hh"
#include "StringUtil.hh"
//#define DEBUG

int GVCFReader::FlushBuffer(bcf1_t *record)
{
    assert(record!=nullptr);
    _depth_buffer.FlushBuffer(record->rid, record->pos - 1);
    int num_flushed = _variant_buffer.FlushBuffer(record);
    FillBuffer();
    return (num_flushed);
}

int GVCFReader::FlushBuffer(int chrom, int pos)
{
    _depth_buffer.FlushBuffer(chrom, pos);
    int num_flushed = _variant_buffer.FlushBuffer(chrom, pos);
    FillBuffer();
    return (num_flushed);
}

int GVCFReader::FlushBuffer()
{
    _depth_buffer.FlushBuffer();
    return (_variant_buffer.FlushBuffer());
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
    FillBuffer();

    // flush variant buffer to get rid of variants overlapping 
    // the interval start
    if(region.find(":")!=std::string::npos)
    {
        string chr;
        int64_t start=0, end = 0;
        stringutil::parsePos(region, chr, start, end);
        if (!region.empty())
        {
            int rid = bcf_hdr_name2id(_bcf_header, chr.c_str());
            FlushBuffer(rid, start);
        }
    }

    //Checking and warning if a few tags are not present. This is how we support legacy GVCFs without crashing.
    if(bcf_hdr_id2int(_bcf_header, BCF_DT_ID, "ADF")==-1)
        std::cerr << "WARNING: " << input_gvcf << " has no FORMAT/ADF tag"<<std::endl;
    if(bcf_hdr_id2int(_bcf_header, BCF_DT_ID, "ADR")==-1)
        std::cerr << "WARNING: " << input_gvcf << " has no FORMAT/ADR tag"<<std::endl;
    if(bcf_hdr_id2int(_bcf_header, BCF_DT_ID, "PL")==-1)
        std::cerr << "WARNING: " << input_gvcf << " has no FORMAT/PL tag" <<std::endl;
    if(bcf_hdr_id2int(_bcf_header, BCF_DT_ID, "MQ")==-1)
        std::cerr << "WARNING: " << input_gvcf << " has no MQ tag" <<std::endl;
}

bool GVCFReader::HasPl()
{
    return(bcf_hdr_id2int(_bcf_header, BCF_DT_ID, "PL")!=-1);
}

bool GVCFReader::HasStrandAd()
{
    return(bcf_hdr_id2int(_bcf_header, BCF_DT_ID, "ADF")!=-1 && bcf_hdr_id2int(_bcf_header, BCF_DT_ID, "ADR")!=-1);
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

size_t GVCFReader::FillBuffer()
{

    if (_variant_buffer.Size() < 2)
    {
        ReadLines(2);
    }
    while (_variant_buffer.Size() > 1 && _variant_buffer.Back()->rid == _variant_buffer.Front()->rid &&
           (_variant_buffer.Back()->pos - _variant_buffer.Front()->pos) < _buffer_size)
    {
        int num_read = ReadLines(1);
        if (num_read == 0)
        {
            break;
        }
    }

    return (_variant_buffer.Size());
}

int GVCFReader::ReadUntil(int rid, int pos)
{
    int num_read = 0;
//    bcf1_t *rec=_variant_buffer.back();
    DepthBlock *db = _depth_buffer.Back();
    while (db == nullptr)
    {
        ReadLines(1);
        db = _depth_buffer.Back();
    }

    while (db != nullptr && bcf_sr_has_line(_bcf_reader, 0) && (db->rid() < rid || (db->rid() == rid && db->end() < pos)))
    {
        if (ReadLines(1) < 1)
        {
            break;
        }
        else
        {
            num_read++;
        }
        db = _depth_buffer.Back();
    }

    return (num_read);
}

int GVCFReader::ReadLines(const unsigned num_lines)
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
            _normaliser->Unarise(_bcf_record, atomised_variants);
            for (auto v = atomised_variants.begin();v!=atomised_variants.end();v++)
            {
                _variant_buffer.PushBack(_bcf_header, *v);
            }
            num_read++;
        }
        int32_t dp;
        //buffer a depth block. FIXME: this should really all be in the DepthBlock constructor.
        if(ggutils::bcf1_get_one_format_int(_bcf_header, _bcf_record, "DP",dp)==1)
        {
            int ploidy = ggutils::get_ploidy(_bcf_header,_bcf_record);
            int start = _bcf_record->pos;
            int32_t dpf, gq, end;
            end = ggutils::get_end_of_gvcf_block(_bcf_header, _bcf_record);
            //If the record has FORMAT/GQ, use that, otherwise take FORMAT/GQX (illumina gvcf quirk).
            int status = ggutils::bcf1_get_one_format_int(_bcf_header,_bcf_record,"GQ",gq);
            if(status!=1)
            {
                float tmp;
                status = ggutils::bcf1_get_one_format_float(_bcf_header, _bcf_record, "GQ", tmp);
                if (status == 1)
                    gq = bcf_float_is_missing(tmp) ? 0 : (int32_t) tmp; //replace missing values with 0
                if (status != 1)
                    status = ggutils::bcf1_get_one_format_int(_bcf_header, _bcf_record, "GQX", gq);
                if (status != 1)
                    ggutils::die("no FORMAT/GQ found");
            }
            gq = gq==bcf_int32_missing ? 0 : gq; //replace missing values with 0
            ggutils::bcf1_get_one_format_int(_bcf_header,_bcf_record,"DPF",dpf);
            _depth_buffer.push_back(DepthBlock(_bcf_record->rid, start, end, dp, dpf, gq, ploidy));
        }
    }
    return (num_read);
}

bcf1_t *GVCFReader::Front()
{
    FillBuffer();
    return (_variant_buffer.Front());
}

bcf1_t *GVCFReader::Pop()
{
    int num_read = FillBuffer();

    if (_variant_buffer.IsEmpty() && num_read == 0)
    {
        return (nullptr);
    }

    return (_variant_buffer.Pop());
}

bool GVCFReader::IsEmpty()
{
    if (_variant_buffer.IsEmpty())
    {
        return (FillBuffer() == 0);
    }
    else
    {
        return (false);
    }
}

bcf_hdr_t *GVCFReader::GetHeader()
{
    return (_bcf_header);
}

//gets dp/dpf/gq (possibly interpolated) for a give interval a<=x<b
void GVCFReader::GetDepth(int rid, int start, int stop, DepthBlock &db)
{
    ReadUntil(rid, stop);
    if (_depth_buffer.Interpolate(rid, start, stop, db) < 0)
    {
        ggutils::die("GVCFReader::get_depth problem with depth buffer");
    }
}

size_t GVCFReader::GetNumVariants()
{
    return (_variant_buffer.Size());
}

size_t GVCFReader::GetNumDepthBlocks()
{
    return (_depth_buffer.size());
}

//gets all variants in interval start<=x<=stop
pair<std::deque<bcf1_t *>::iterator,std::deque<bcf1_t *>::iterator> GVCFReader::GetAllVariantsInInterval(int chrom,
                                                                                                         int stop)
{
    return(_variant_buffer.GetAllVariantsInInterval(chrom, stop));
}

pair<std::deque<bcf1_t *>::iterator,std::deque<bcf1_t *>::iterator>  GVCFReader::GetAllVariantsUpTo(bcf1_t *record)
{
    return(_variant_buffer.GetAllVariantsUpTo(record));
}

