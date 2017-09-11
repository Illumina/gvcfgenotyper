#include "GVCFReader.hpp"
//#define DEBUG

static void remove_hdr_lines(bcf_hdr_t *hdr, int type)
{
    int i = 0, nrm = 0;
    while ( i<hdr->nhrec )
    {
        if ( hdr->hrec[i]->type!=type ) 
	{ 
	    i++; 
	    continue; 
	}
        bcf_hrec_t *hrec = hdr->hrec[i];
        if ( type==BCF_HL_FMT )
        {
            // everything except FORMAT/GT
            int id = bcf_hrec_find_key(hrec, "ID");
            if ( id>=0 && !strcmp(hrec->vals[id],"GT") ) 
	    {
		i++; continue;
	    }
        }
        nrm++;
        hdr->nhrec--;
        if ( i < hdr->nhrec )
            memmove(&hdr->hrec[i],&hdr->hrec[i+1],(hdr->nhrec-i)*sizeof(bcf_hrec_t*));
        bcf_hrec_destroy(hrec);
    }
    if ( nrm ) bcf_hdr_sync(hdr);
}

void remove_info(bcf1_t *line)
{
    // remove all INFO fields
    if ( !(line->unpacked & BCF_UN_INFO) ) bcf_unpack(line, BCF_UN_INFO);

    int i;
    for (i=0; i<line->n_info; i++)
    {
        bcf_info_t *inf = &line->d.info[i];
        if ( inf->vptr_free )
        {
            free(inf->vptr - inf->vptr_off);
            inf->vptr_free = 0;
        }
        line->d.shared_dirty |= BCF1_DIRTY_INF;
        inf->vptr = NULL;
    }
    line->n_info=0;
}

int GVCFReader::flush_buffer(const bcf1_t *record)
{
    return(_variant_buffer.flush_buffer(record));
}

int GVCFReader::flush_buffer(int chrom,int pos)
{
    return(_variant_buffer.flush_buffer(chrom,pos));
}  


int GVCFReader::flush_buffer()
{
    return(_variant_buffer.flush_buffer());
}  

GVCFReader::GVCFReader(const std::string & input_gvcf,const std::string & reference_genome_fasta,int buffer_size)  
{
    _bcf_record=NULL;
    _bcf_reader =  bcf_sr_init() ; 
    if(!(bcf_sr_add_reader (_bcf_reader, input_gvcf.c_str())))
    {
	die("problem opening input");
    }
    if(buffer_size<2)
    {
	die("GVCFReader needs buffer size of at least 2");
    }
    _buffer_size=buffer_size;
    _bcf_header= _bcf_reader->readers[0].header;    
    _bcf_record=NULL;
    //this is a hack to fix gvcfs where AD is set to Number=. (VCF4.1 does not technically allow Number=R)
    bcf_hdr_remove(_bcf_header,BCF_HL_FMT,"AD");
    assert(  bcf_hdr_append(_bcf_header,"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed. For indels this value only includes reads which confidently support each allele (posterior prob 0.999 or higher that read contains indicated allele vs all other intersecting indel alleles)\">") == 0);

    //this is a hack to fix gvcfs where GQ is labelled as float (VCF spec says it should be integer)
    bcf_hdr_remove(_bcf_header,BCF_HL_FMT,"GQ");
    assert(  bcf_hdr_append(_bcf_header,"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">") == 0);
    _normaliser = new Normaliser(reference_genome_fasta,_bcf_header);
    fill_buffer(buffer_size);
}

GVCFReader::~GVCFReader()
{
    if(_bcf_record!=NULL)
    {
	bcf_destroy1(_bcf_record);
    }
    if ( _bcf_reader->errnum )
    {
	error("Error: %s\n", bcf_sr_strerror(_bcf_reader->errnum));
    }
//    bcf_sr_destroy(_bcf_reader); //this is causing an invalid free. i am not sure why!
    delete _normaliser;    
}

int GVCFReader::fill_buffer(int num_lines)
{
    return(read_lines(num_lines - _variant_buffer.size()));
}

int GVCFReader::read_lines(int num_lines)
{
    if(num_lines<=0)
    {
	return(0);
    }

    int num_read=0;
   
    while(num_read<num_lines && bcf_sr_next_line(_bcf_reader))
    {
	_bcf_record =  bcf_sr_get_line(_bcf_reader, 0);
   	if(_bcf_record->n_allele>1)//is this line a variant?
	{
	    int32_t pass = bcf_has_filter(_bcf_header, _bcf_record, (char *)".");
	    bcf_update_format_int32(_bcf_header,_bcf_record,"FT",&pass,1);
	    bcf_update_filter(_bcf_header,_bcf_record,NULL,0);
	    bcf_update_id(_bcf_header,_bcf_record,NULL);
	    remove_info(_bcf_record);
	    vector<bcf1_t *> atomised_variants = _normaliser->atomise(_bcf_record);
	    for(size_t i=0;i<atomised_variants.size();i++)
	    {
//		print_variant(_bcf_header,atomised_variants[i]);//debug
		_variant_buffer.push_back(atomised_variants[i]);
//		bcf_destroy1(atomised_variants[i]);// VariantBuffer handles memory.
	    }
	    num_read++;
	}

	//buffer a depth block
	int num_format_values=0;
	int32_t *value_pointer=NULL;
//if DP is present, this is either a snp or a homref block and we want to store it in depth buffer;
	if(bcf_get_format_int32(_bcf_header, _bcf_record, "DP", &value_pointer,&num_format_values)==1)
	{
	    int start=_bcf_record->pos;
	    int32_t dp,dpf,gq,end;
	    dp = *value_pointer;
	    end = get_end_of_gvcf_block(_bcf_header,_bcf_record);
	    //if it is a SNP use GQ else use GQX (this is a illumina GVCF quirk)
	    if(_bcf_record->n_allele>1)
	    {
		bcf_get_format_int32(_bcf_header, _bcf_record, "GQ", &value_pointer , &num_format_values);
	    }	    
	    else
	    {
		bcf_get_format_int32(_bcf_header, _bcf_record, "GQX", &value_pointer , &num_format_values);
	    }
	    gq = *value_pointer;
	    bcf_get_format_int32(_bcf_header, _bcf_record, "DPF", &value_pointer , &num_format_values);
	    dpf = *value_pointer;
	    _depth_buffer.push_back(DepthBlock(_bcf_record->rid,start,end,dp,dpf,gq));
	    free(value_pointer);
	}

    }

    return(num_read);
}

bcf1_t *GVCFReader::front()
{
    if(_variant_buffer.empty())
    {
	read_lines(_buffer_size);
    }
    return(_variant_buffer.front());
}

bcf1_t *GVCFReader::pop()
{
    int num_read = 0;
    if(_variant_buffer.size() < _buffer_size/2)
    {
	read_lines(_buffer_size-_variant_buffer.size());
    }

    if(_variant_buffer.empty() && num_read==0)
    {
	return(NULL);
    }

    return(_variant_buffer.pop());
}

bool GVCFReader::empty()
{
    if(_variant_buffer.empty())
    {
	return(read_lines(_buffer_size)==0);
    }
    else
    {
	return(false);
    }
}

const bcf_hdr_t *GVCFReader::get_header()
{
    return(_bcf_header);
}

//gets dp/dpf/gq (possibly interpolated) for a give interval a<=x<b
void GVCFReader::get_depth(int rid,int start,int stop,DepthBlock & db)
{
    _depth_buffer.interpolate(rid,start,stop,db);
}
