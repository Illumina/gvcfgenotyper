#include "GVCFReader.hpp"
//#define DEBUG

using namespace std;

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


Normaliser::Normaliser(const string & ref_fname,bcf_hdr_t *hdr)
{
    _hdr=bcf_hdr_dup(hdr);
    _norm_args=init_vcfnorm(_hdr,ref_fname.c_str());
}

Normaliser::~Normaliser()
{
    bcf_hdr_destroy(_hdr);
    destroy_data(_norm_args);
    free(_norm_args);
}

//1. split multi-allelics
//2. normalise (left-align + trim)
//3. decompose MNPs into SNPs
vector<bcf1_t *>  Normaliser::atomise(bcf1_t *bcf_record_to_canonicalise)
{
    assert(bcf_record_to_canonicalise->n_allele>1);
    vector<bcf1_t *> atomised_variants;
    bcf1_t **split_records=&bcf_record_to_canonicalise;
    int num_split_records=1;
    if(bcf_record_to_canonicalise->n_allele>2)
    {//split multi-allelics (using vcfnorm.c from bcftools1.3
	split_multiallelic_to_biallelics(_norm_args,bcf_record_to_canonicalise);
	split_records=_norm_args->tmp_lines;
	num_split_records=_norm_args->ntmp_lines;
    }
	    
    for(int rec_index=0;rec_index<num_split_records;rec_index++)
    {
	bcf1_t *rec=split_records[rec_index];
//	cerr << "atomise: "<<rec->pos+1<<":"<<rec->d.allele[0]<<":"<<rec->d.allele[1]<<endl;
	if(strlen(rec->d.allele[0])==1 && strlen(rec->d.allele[1])==1)//is a snp. do nothing
	{
	    atomised_variants.push_back(bcf_dup(rec));
	}
	else
	{
	    if(realign(_norm_args,rec) == ERR_REF_MISMATCH)
	    {
		die("vcf record did not match the reference");
	    }
	    char *ref=rec->d.allele[0];
	    char *alt=rec->d.allele[1];
	    int ref_len = strlen(ref);
	    int alt_len = strlen(alt);

	    if(ref_len>1 && ref_len==alt_len) //is MNP
	    {
		char alleles[4] = "X,X";
		for(int i=0;i<ref_len;i++) 
		{
		    if(ref[i]!=alt[i]) 
		    {//new SNP
			bcf1_t *new_var = bcf_dup(rec);
			bcf_unpack(new_var, BCF_UN_ALL);
			alleles[0]=ref[i];
			alleles[2]=alt[i];
			new_var->pos+=i;
			bcf_update_alleles_str(_hdr, new_var, alleles);	
			atomised_variants.push_back(new_var);
//			cerr << "new_var: "<<new_var->pos+1<<":"<<new_var->d.allele[0]<<":"<<new_var->d.allele[1]<<endl;
		    }
		}
	    }
	    // else if((ref_len!=alt_len) && (ref_len!=1) && (alt_len>1)) //complex substitution
	    // {
	    //     vt_aggressive_decompose(rec,_hdr,atomised_variants);
	    // }
	    else //variant already is atomic
	    {
		bcf1_t *new_var = bcf_dup(rec);
		atomised_variants.push_back(new_var);
	    }  
	}
    }
    return(atomised_variants);
}

VariantBuffer::VariantBuffer() 
{
    _num_duplicated_records=0;
}

VariantBuffer::~VariantBuffer() 
{
    cerr << "Dropped "<<_num_duplicated_records <<" duplicated variants after normalization."<<endl;
    flush_buffer();
}

bool VariantBuffer::has_variant(bcf1_t *v) 
{
    int i = _buffer.size()-1;
    while(i>=0 && _buffer[i]->pos >= v->pos)
    {
	if(v==_buffer[i])
	{
	    return(true);
	}
	i--;
    }
    return(false);
}

int VariantBuffer::push_back(bcf1_t *v) 
{
    bcf_unpack(v, BCF_UN_ALL);
    if(has_variant(v))
    {
	_num_duplicated_records++;
	return(0);
    }

    _buffer.push_back(bcf_dup(v));
    int i = _buffer.size()-1;
    while(i>0 && _buffer[i]->pos < _buffer[i-1]->pos) 
    {
	bcf1_t *tmp=_buffer[i-1];
	_buffer[i-1]=_buffer[i];
	_buffer[i]=tmp;
	i--;
    }
    return(1);
}

int GVCFReader::flush_buffer(int chrom,int pos)
{
    return(_variant_buffer.flush_buffer(chrom,pos));
}  


int GVCFReader::flush_buffer()
{
    return(_variant_buffer.flush_buffer());
}  


int VariantBuffer::flush_buffer(int chrom,int pos)
{
    int num_flushed=0;
    while(_buffer.size()>0 &&  _buffer.front()->pos <= pos && _buffer.front()->rid <= chrom)
    {
	bcf1_t *rec = _buffer.front();	    
	_buffer.pop_front();
	num_flushed++;
    }    
    return(num_flushed);
}  

int VariantBuffer::flush_buffer()
{
    if(_buffer.size()>0) 
    {
	int ret = flush_buffer(_buffer.back()->rid,_buffer.back()->pos);
	return(ret);
    }
    else
    {
	return(0);
    }
}

GVCFReader::GVCFReader(const std::string & input_gvcf,const std::string & reference_genome_fasta,int buffer_size)  
{
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
    delete _normaliser;
//    bcf_sr_destroy(_bcf_reader);
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
	    int32_t pass = bcf_has_filter(_bcf_header, _bcf_record, ".");
	    bcf_update_format_int32(_bcf_header,_bcf_record,"FT",&pass,1);
	    bcf_update_filter(_bcf_header,_bcf_record,NULL,0);
	    bcf_update_id(_bcf_header,_bcf_record,NULL);
	    remove_info(_bcf_record);
	    vector<bcf1_t *> atomised_variants = _normaliser->atomise(_bcf_record);
	    for(size_t i=0;i<atomised_variants.size();i++)
	    {
		_variant_buffer.push_back(atomised_variants[i]);
		bcf_destroy1(atomised_variants[i]);
	    }
	    num_read++;
	}
    }

    return(num_read);
}

bcf1_t *GVCFReader::front()
{
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

size_t VariantBuffer::size()
{
    return(_buffer.size());
}

bool VariantBuffer::empty()
{
    return(_buffer.empty());
}

bcf1_t *VariantBuffer::front()
{
    if(_buffer.empty())
    {
	return(NULL);
    }
    else
    {	
	bcf1_t *ret = _buffer.front();
	bcf_unpack(ret, BCF_UN_ALL);
	return(ret);
    }
}


bcf1_t *VariantBuffer::pop()
{
    if(_buffer.empty())
    {
	return(NULL);
    }
    else
    {
	bcf1_t *ret = _buffer.front();
	bcf_unpack(ret, BCF_UN_ALL);
	_buffer.pop_front();
	return(ret);
    }
}

const bcf_hdr_t *GVCFReader::getHeader()
{
    return(_bcf_header);
}
