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


Normaliser::Normaliser(const string & ref_fname)
{
    _norm_args  = (args_t*) calloc(1,sizeof(args_t));
  _norm_args->files   = NULL;
  _norm_args->output_fname = NULL;
  _norm_args->output_type = FT_VCF;
  _norm_args->aln_win = 100;
  _norm_args->buf_win = 1000;
  _norm_args->mrows_collapse = COLLAPSE_BOTH;
  _norm_args->mrows_op = MROWS_SPLIT;
  _norm_args->hdr = NULL;//hdr;
  _norm_args->do_indels = 1;
  _norm_args->ref_fname = (char *)ref_fname.c_str();
  init_data(_norm_args);
}


//1. decompose MNPs/complex substitutions 
//2. normalises the output using bcftools norm algorithm
vector<bcf1_t *>  Normaliser::atomise(bcf1_t *rec,bcf_hdr_t *hdr)
{
    assert(rec->n_allele>1);
    vector<bcf1_t *> atomised_variants;
    bcf1_t **split_records=&rec;
    int num_split_records=1;
    if(rec->n_allele>2)
    {//split multi-allelics (using vcfnorm.c from bcftools1.3
	split_multiallelic_to_biallelics(_norm_args,rec);
	split_records=_norm_args->tmp_lines;
	num_split_records=_norm_args->ntmp_lines;
    }
	    
    for(int i=0;i<num_split_records;i++)
    {
	bcf1_t *rec=split_records[i];
	char *ref=rec->d.allele[0];
	char *alt=rec->d.allele[1];
	int ref_len = strlen(ref);
	int alt_len = strlen(alt);
	if(ref_len==1 && alt_len==1)//it is a SNP. no more work required.
	{
	    atomised_variants.push_back(bcf_dup(rec));
	}
	else
	{
	    if(realign(_norm_args,rec) == ERR_REF_MISMATCH)
	    {
		die("vcf record did not match the reference");
	    }
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
			bcf_update_alleles_str(hdr, new_var, alleles);	
			atomised_variants.push_back(new_var);
		    }
		}
	    }
	    // else if((ref_len!=alt_len) && (ref_len!=1) && (alt_len>1)) //complex substitution
	    // {
	    //     vt_aggressive_decompose(rec,hdr,atomised_variants);
	    // }
	    else //variant already is atomic
	    {
		bcf1_t *new_var = bcf_dup(rec);
		bcf_unpack(new_var, BCF_UN_ALL);
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

GVCFReader::GVCFReader(const std::string & input_gvcf,const std::string & reference_genome_fasta)  
{
    _bcf_reader =  bcf_sr_init() ; 
    _bcf_header= _bcf_reader->readers[0].header;
    if(!(bcf_sr_add_reader (_bcf_reader, input_gvcf.c_str())))
    {
	die("problem opening input");
    }
    
    //this is a hack to fix gvcfs where AD is set to Number=. (VCF4.1 does not technically allow Number=R)
    bcf_hdr_remove(_bcf_header,BCF_HL_FMT,"AD");
    assert(  bcf_hdr_append(_bcf_header,"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed. For indels this value only includes reads which confidently support each allele (posterior prob 0.999 or higher that read contains indicated allele vs all other intersecting indel alleles)\">") == 0);

    //this is a hack to fix gvcfs where GQ is labelled as float (VCF spec says it should be integer)
    bcf_hdr_remove(_bcf_header,BCF_HL_FMT,"GQ");
    assert(  bcf_hdr_append(_bcf_header,"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">") == 0);
    _normaliser = new Normaliser(reference_genome_fasta);
}

GVCFReader::~GVCFReader()
{
    bcf_destroy1(_bcf_record);
}

int GVCFReader::read_lines(int num_lines)
{
    int num_read=0;
   
    while(bcf_sr_next_line(_bcf_reader) && num_read<num_lines)
    {
	_bcf_record =  bcf_sr_get_line(_bcf_reader, 0);
   	if(_bcf_record->n_allele>1)//is this line a variant?
	{
	    int32_t pass = bcf_has_filter(_bcf_header, _bcf_record, ".");
	    bcf_update_format_int32(_bcf_header,_bcf_record,"FT",&pass,1);
	    bcf_update_filter(_bcf_header,_bcf_record,NULL,0);
	    bcf_update_id(_bcf_header,_bcf_record,NULL);
	    remove_info(_bcf_record);
	    vector<bcf1_t *> atomised_variants = _normaliser->atomise(_bcf_record,_bcf_header);
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

