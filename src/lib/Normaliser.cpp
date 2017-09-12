#include "GVCFReader.hpp"
#include "needle.hpp"
#include "vt_utils.hpp"
#include "hts_utils.hpp"

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
//	std::cerr << "atomise: "<<rec->pos+1<<":"<<rec->d.allele[0]<<":"<<rec->d.allele[1]<<std::endl;
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
//			std::cerr << "new_var: "<<new_var->pos+1<<":"<<new_var->d.allele[0]<<":"<<new_var->d.allele[1]<<std::endl;
		    }
		}
	    }
	    else if((ref_len!=alt_len) && (ref_len!=1) && (alt_len>1)) //complex substitution
	    {
	        vt_aggressive_decompose(rec,_hdr,atomised_variants);
	    }
	    else //variant already is atomic
	    {
		bcf1_t *new_var = bcf_dup(rec);
		atomised_variants.push_back(new_var);
	    }  
	}
    }
    return(atomised_variants);
}
