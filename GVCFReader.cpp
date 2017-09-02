#include "agg_ingest1.h"
//#define DEBUG

//just a struct to count some things
struct Counts
{
    int mnp,complex;
};

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

char *find_format(char *ptr,char *FORMAT) 
{
    char *fmt_ptr = strstr(ptr,FORMAT);
    if(fmt_ptr!=NULL) 
    {
	int idx=0;
	for(int i=0;i<(fmt_ptr-ptr);i++)
	    if(ptr[i]==':')
		idx++;	
	ks_tokaux_t aux;
	ptr = kstrtok(ptr,"\t",&aux);
	ptr = kstrtok(NULL,NULL,&aux);//sample column
	ptr = kstrtok(ptr,":",&aux);
	for(int i=0;i<idx;i++) 
	    ptr = kstrtok(NULL,NULL,&aux);

	return(ptr);
    }
    else
	return(NULL);
}


//vt triple structure
struct Triple
{
    int pos_ref;
    int pos_alt;
    int len_ref;
    int len_alt;

    Triple() : pos_ref(0), pos_alt(0), len_ref(0), len_alt(0)
    {}

    Triple(int pos_ref, int pos_alt, int len_ref, int len_alt) : pos_ref(pos_ref), pos_alt(pos_alt), len_ref(len_ref), len_alt(len_alt)
    {}
};

//this is a (slightly) modfified version of vt's aggressive decompose 
//see https://github.com/atks/vt
int vt_aggressive_decompose(bcf1_t *v,bcf_hdr_t *hdr,vector<bcf1_t *> & buf)
{


    char** allele = bcf_get_allele(v);


    int new_no_variants=0;
    kstring_t new_alleles= {0,0,0};
    // Use alignment for decomposition of substitutions where REF
    // and ALT have different lengths and the variant is not an
    // insertion or deletion.
    
    // Perform alignment of REF[1:] and ALT[1:]
    NeedlemanWunsch nw(true);
    nw.align(allele[0] + 1, allele[1] + 1);
    nw.trace_path();
    // Force-align first characters
    if (allele[0][0] == allele[1][0])
	nw.trace.insert(nw.trace.begin(), NeedlemanWunsch::CIGAR_M);
    else
	nw.trace.insert(nw.trace.begin(), NeedlemanWunsch::CIGAR_X);
    nw.read--;
    nw.ref--;

    // Break apart alignment
    std::vector<Triple> chunks;
    bool hasError = false;
    int pos_ref = 0, pos_alt = 0, k = 0;
    Triple nextChunk(pos_ref, pos_alt, 0, 0);
    while (pos_ref <= nw.len_ref || pos_alt <= nw.len_read)
    {
	switch ((int32_t)nw.trace.at(k++))
	{
	case NeedlemanWunsch::CIGAR_M:
	    if (hasError)
		chunks.push_back(nextChunk);
	    nextChunk = Triple(pos_ref++, pos_alt++, 1, 1);
	    hasError = false;
	    break;
	case NeedlemanWunsch::CIGAR_X:
	    if (hasError)
		chunks.push_back(nextChunk);
	    nextChunk = Triple(pos_ref++, pos_alt++, 1, 1);
	    hasError = true;
	    break;
	case NeedlemanWunsch::CIGAR_D:
	    nextChunk.len_ref++;
	    pos_ref++;
	    hasError = true;
	    break;
	case NeedlemanWunsch::CIGAR_I:
	    nextChunk.len_alt++;
	    pos_alt++;
	    hasError = true;
	    break;
	}
    }
    if (hasError)
	chunks.push_back(nextChunk);


    int32_t pos1 = bcf_get_pos1(v);
    char* ref = strdup(v->d.allele[0]);
    char* alt = strdup(v->d.allele[1]);

    // old_alleles.l = 0;
    // bcf_variant2string(hdr, v, &old_alleles);

    for (size_t i=0; i<chunks.size(); ++i)
    {
	bcf1_t *nv = bcf_dup(v);
	bcf_unpack(nv, BCF_UN_ALL);
	bcf_set_pos1(nv, pos1+chunks[i].pos_ref);
	std::vector<int32_t> start_pos_of_phased_block;

                        
	new_alleles.l=0;
	for (int j=chunks[i].pos_ref; j<chunks[i].pos_ref+chunks[i].len_ref; ++j)
	    kputc(ref[j], &new_alleles);
	kputc(',', &new_alleles);
	for (int j=chunks[i].pos_alt; j<chunks[i].pos_alt+chunks[i].len_alt; ++j)
	    kputc(alt[j], &new_alleles);

	bcf_update_alleles_str(hdr, nv, new_alleles.s);
//	bcf_update_info_string(hdr, nv, "OLD_CLUMPED", old_alleles.s);
                    
	buf.push_back(nv);
	kputc('\0', &new_alleles);

	++new_no_variants;
    }
    if(new_alleles.l)
    {
	free(new_alleles.s);
    }

    free(ref);
    free(alt);
    
    return(new_no_variants);
}


//1. decompose MNPs/complex substitutions 
//2. normalises the output using bcftools norm algorithm
vector<bcf1_t *> atomise(bcf1_t *rec,bcf_hdr_t *hdr,Counts & counts)
{
    assert(rec->n_allele  == 2);
    char *ref=rec->d.allele[0];
    char *alt=rec->d.allele[1];
    int ref_len = strlen(ref);
    int alt_len = strlen(alt);
    vector<bcf1_t *> ret;
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
		ret.push_back(new_var);
	    }
	    counts.mnp++;	    
	}
    }
    else if((ref_len!=alt_len) && (ref_len!=1) && (alt_len>1)) //complex substitution
    {
	vt_aggressive_decompose(rec,hdr,ret);
	counts.complex++;
    }
    else //variant already is atomic
    {
	bcf1_t *new_var = bcf_dup(rec);
	bcf_unpack(new_var, BCF_UN_ALL);
	ret.push_back(new_var);
    }  
    return(ret);
}

int ingest1(const char *input,const char *output,char *ref,bool exit_on_mismatch=true) 
{
    cerr << "Input: " << input << "\tOutput: "<<output<<endl;
    Counts counts;
    counts.mnp=0;
    counts.complex=0;
    gzFile fp = gzopen(input, "r");
    VarBuffer vbuf(1000);
    int prev_rid = -1;
    if(fp==NULL) {
	fprintf(stderr,"problem opening %s\n",input);
	exit(1);
    }

    char *out_fname = (char *)malloc(strlen(output)+5);
    strcpy(out_fname,output);
    strcat(out_fname,".tmp");
    if(fileexists(out_fname)) {
	fprintf(stderr,"%s file already exists. will not overwrite\n",out_fname);
	exit(1);
    }
    printf("depth: %s\n",out_fname);
    gzFile depth_fp = gzopen(out_fname, "wb1");
    strcpy(out_fname,output);
    strcat(out_fname,".bcf");
    if(fileexists(out_fname)) {
	fprintf(stderr,"%s file already exists. will not overwrite\n",out_fname);
	exit(1);
    }
    printf("variants: %s\n",out_fname);
    htsFile *variant_fp=hts_open(out_fname,"wb1");
    if(variant_fp==NULL) {
	fprintf(stderr,"problem opening %s\n",input);
	exit(1);    
    }

    ks = ks_init(fp);
    htsFile *hfp=hts_open(input, "r");
    bcf_hdr_t *hdr_in =  bcf_hdr_read(hfp);
    if(bcf_hdr_id2int(hdr_in, BCF_DT_ID, "GQX")==-1)
	die("FORMAT/GQX is not present");
    if(bcf_hdr_id2int(hdr_in, BCF_DT_ID, "DP")==-1)
	die("FORMAT/DP is not present");
    if(bcf_hdr_id2int(hdr_in, BCF_DT_ID,"BLOCKAVG_min30p3a")==-1)
	die("INFO/BLOCKAVG_min30p3a is not present");

    hts_close(hfp);
    //this is a hack to fix gvcfs where AD is incorrectly defined in the header. (vcf4.2 does not technically allow Number=R)
    bcf_hdr_remove(hdr_in,BCF_HL_FMT,"AD");
    assert(  bcf_hdr_append(hdr_in,"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed. For indels this value only includes reads which confidently support each allele (posterior prob 0.999 or higher that read contains indicated allele vs all other intersecting indel alleles)\">") == 0);

    //this is a hack to fix broken gvcfs where GQ is incorrectly labelled as float (v4.3 spec says it should be integer)
    bcf_hdr_remove(hdr_in,BCF_HL_FMT,"GQ");
    assert(  bcf_hdr_append(hdr_in,"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">") == 0);


    //  bcf_hdr_t  *hdr_out=hdr_in;
    bcf_hdr_t *hdr_out =  bcf_hdr_dup(hdr_in);
    remove_hdr_lines(hdr_out,BCF_HL_INFO);
    remove_hdr_lines(hdr_out,BCF_HL_FLT);
    bcf_hdr_sync(hdr_out);

    //here we add FORMAT/PF. which is the pass filter flag for alts.
    assert(  bcf_hdr_append(hdr_out,"##FORMAT=<ID=PF,Number=A,Type=Integer,Description=\"variant was PASS filter in original sample gvcf\">") == 0);

    args_t *norm_args = init_vcfnorm(hdr_out,ref);
    norm_args->check_ref |= CHECK_REF_WARN;
    bcf_hdr_write(variant_fp, hdr_out);

    int buf[5];

    bcf_srs_t *sr =  bcf_sr_init() ; ///htslib synced reader.
    if(!(bcf_sr_add_reader (sr, input )))
    {
	die("problem opening input");
    }
    bcf1_t *line;
    int dp,gq,end;
    int *dp_ptr = &dp;
    int *gq_ptr = &gq;
    int *end_ptr = &end;
    int nval=1;
    while(bcf_sr_next_line (sr))
    {
	line =  bcf_sr_get_line(sr, 0);
	bcf_unpack(line,BCF_UN_ALL);
	if(bcf_get_format_int32(sr->readers[0].header, line, "DP", &dp_ptr , &nval) == 1)
	{
	    buf[0] = line->rid;
	    buf[1] = line->pos;
	    if(bcf_get_format_int32(sr->readers[0].header, line, "END", &end_ptr , &nval) == 1)
	    {
		buf[2] = end - 1;
	    }
	    else
	    {
		int alt_len = 1;
		if(line->n_allele>1)
		{
		    alt_len = strlen(line->d.allele[1]);
		}
		buf[2] = buf[1]+alt_len-1;
	    }
	    buf[3] = dp;
	    if(bcf_get_format_int32(sr->readers[0].header, line, "GQ", &gq_ptr , &nval) == 1)
	    {
		buf[4] = gq;
	    }
	    else
	    {
		assert(bcf_get_format_int32(sr->readers[0].header, line, "GQX", &gq_ptr , &nval)==1);
		buf[4]=gq/10;
		buf[4]*=10;
	    }

#ifdef DEBUG
		fprintf(stderr,"%d\t%d\t%d\t%d\t%d\n",buf[0],buf[1],buf[2],buf[3],buf[4]);
#endif 
		if(gzwrite(depth_fp,buf,5*sizeof(int))!=(5*sizeof(int)))
		{
		    die("ERROR: problem writing "+(string)out_fname+".tmp");		    
		}
	}

	if(line->n_allele>1)
	{//was this a variant? if so write it out to the bcf
	    norm_args->ntotal++;
	    //	cerr<<line->rid<<":"<<line->pos<<endl;
	    if(prev_rid!=line->rid)
	    {
		vbuf.flush(variant_fp,hdr_out);
	    }		 
	    else
	    {
		vbuf.flush(line->pos,variant_fp,hdr_out);
	    }		    
	    prev_rid=line->rid;
	    int32_t pass = bcf_has_filter(hdr_in, line, ".");
	    bcf_update_format_int32(hdr_out,line,"PF",&pass,1);
	    bcf_update_filter(hdr_out,line,NULL,0);
	    bcf_update_id(hdr_out,line,NULL);
	    bcf1_t **split_records=&line;
	    int num_split_records=1;
	    if(line->n_allele>2)
	    {//split multi-allelics (using vcfnorm.c from bcftools1.3
		norm_args->nsplit++;
		split_multiallelic_to_biallelics(norm_args,line);
		split_records=norm_args->tmp_lines;
		num_split_records=norm_args->ntmp_lines;
	    }
	    
	    for(int i=0;i<num_split_records;i++)
	    {
		remove_info(split_records[i]);
		vector<bcf1_t *> atomised_variants = atomise(split_records[i],hdr_out,counts);
		for(size_t j=0;j<atomised_variants.size();j++)
		{
		    if(realign(norm_args,atomised_variants[j]) == ERR_REF_MISMATCH)
		    {
			if(exit_on_mismatch)
			{
				die("vcf did not match the reference");
			}
			else
			{
			    norm_args->nskipped++;
			}			    
		    }
		    vbuf.push_back(atomised_variants[j]);
//			bcf_destroy1(atomised_variants[j]);
		}		    
	    }
	    vbuf.flush(line->pos,variant_fp,hdr_out);
	}
    }

    vbuf.flush(variant_fp,hdr_out);
    bcf_hdr_destroy(hdr_in);
    bcf_hdr_destroy(hdr_out);
    bcf_destroy1(line);
    ks_destroy(ks);
    gzclose(fp);
    gzclose(depth_fp);  

    hts_close(variant_fp);

    fprintf(stderr,"Variant lines   total/split/realigned/skipped:\t%d/%d/%d/%d\n", norm_args->ntotal,norm_args->nsplit,norm_args->nchanged,norm_args->nskipped);
    fprintf(stderr,"Decomposed %d MNPs\n", counts.mnp);
    fprintf(stderr,"Decomposed %d complex substitutions\n", counts.complex);
    
    destroy_data(norm_args);
    free(norm_args);


    fprintf(stderr,"Indexing %s\n",out_fname);
    bcf_index_build(out_fname, BCF_LIDX_SHIFT);
    free(out_fname);
    return 0;
}

VariantBuffer::VariantBuffer() 
{
    _num_duplicated_records=0;
}

VariantBuffer::~VariantBuffer() 
{
    cerr << "Dropped "<<_num_duplicated_records <<" duplicated variants after normalization."<<endl;
    flush();
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


int Variant::flush_buffer(int chrom,int pos)
{
    while(_buffer.size()>0 &&  _buffer.front()->pos <= pos && _buffer.front()->rid <= rid)
    {
	bcf1_t *rec = _buffer.front();	    
	_buffer.pop_front();
    }    
    return(n);
}  

int flush_buffer()
{
    if(_buffer.size()>0) 
    {
	int ret = flush(_buffer.back()->pos);
	return(ret);
    }
    else
    {
	return(0);
    }
}
  
