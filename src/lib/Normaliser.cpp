#include <htslib/vcf.h>
#include "GVCFReader.hpp"

Normaliser::Normaliser(const string &ref_fname, bcf_hdr_t *hdr)
{
    _hdr = hdr;
    _norm_args = init_vcfnorm(_hdr, ref_fname.c_str());
}

Normaliser::~Normaliser()
{
    destroy_data(_norm_args);
    free(_norm_args);
}

int mnp_split(bcf1_t *record_to_split,bcf_hdr_t *header,vector<bcf1_t *> & output)
{
    char *ref = record_to_split->d.allele[0];
    char *alt = record_to_split->d.allele[1];
    size_t ref_len = strlen(ref);
    size_t alt_len = strlen(alt);
    int num_new_snps=0;
    if (ref_len > 1 && ref_len == alt_len) //is MNP
    {
        char alleles[4] = "X,X";
        for (int i = 0; i < ref_len; i++)
        {
            if (ref[i] != alt[i])
            {//new SNP
                bcf1_t *new_var = bcf_dup(record_to_split);
                bcf_unpack(new_var, BCF_UN_ALL);
                alleles[0] = ref[i];
                alleles[2] = alt[i];
                new_var->pos += i;
                bcf_update_alleles_str(header, new_var, alleles);
                output.push_back(new_var);
                num_new_snps++;
            }
        }
        return(num_new_snps);
    }
    else
    {
        output.push_back(bcf_dup(record_to_split));
        return(1);
    }
}

vector<bcf1_t *> Normaliser::unarise(bcf1_t *bcf_record_to_marginalise)
{
    bcf_unpack(bcf_record_to_canonicalise, BCF_UN_ALL);
    assert(bcf_record_to_canonicalise->n_allele > 1);
    vector<bcf1_t *> atomised_variants;

    //this chunk of codes reads the relevant FORMAT fields (PL,GQ,DP,DPF,AD) that will be propagated into the new record.
    int32_t *old_ad,*old_gt,*old_gq,*old_dpf,*old_pl,*old_dp;
    int num_ad=0,num_gt=0,num_gq=0,num_dpf=0,num_pl=0,num_dp=0;

    int ploidy = bcf_get_genotypes(_hdr,bcf_record_to_marginalise,&old_gt,&num_gt);
    assert(ploidy>0);
    assert(bcf_get_format_int32(_hdr,bcf_record_to_marginalise,"AD",&old_ad,&num_ad)>0);
    assert(bcf_get_format_int32(_hdr,bcf_record_to_marginalise,"DPF",&old_dpf,&num_dpf)>0);
    assert(bcf_get_format_int32(_hdr,bcf_record_to_marginalise,"DP",&old_dp,&num_dp)>0);
    if(bcf_get_format_int32(_hdr,bcf_record_to_marginalise,"GQ",&old_gq,&num_gq)==-2)//catches cases where GQ is set as a float (should be int according to vcf4.3)
    {
        float *tmp_gq;
        assert(bcf_get_format_float(_hdr,bcf_record_to_marginalise,"GQ",&tmp_gq,&num_gq)>0);
        old_gq = new int32_t[1];
        old_gq[0] = (int32_t)tmp_gq[0];
        delete tmp_gq;
    }

    //this block deals with bi-allelic variants. no marginalisation across other alternate alleles is required
    if(bcf_record_to_marginalise->n_allele == 2)
    {
        bcf1_t *new_record = bcf_init1();
        new_record->rid = bcf_record_to_marginalise->rid;
        new_record->pos = bcf_record_to_marginalise->pos;
        bcf_update_alleles(_hdr,new_record,bcf_record_to_marginalise->d.allele,bcf_record_to_marginalise->n_allele);
        bcf_update_genotypes(_hdr,new_record,old_gt,num_gt);
        bcf_update_format_int32(_hdr,new_record,"GQ",old_gq,num_gq);
        bcf_update_format_int32(_hdr,new_record,"AD",old_ad,num_ad);
        bcf_update_format_int32(_hdr,new_record,"DP",old_dp,num_dp);
        bcf_update_format_int32(_hdr,new_record,"DPF",old_dpf,num_dpf);
        bcf_update_format_int32(_hdr,new_record,"PL",old_pl,num_pl);

        if(strlen(bcf_record_to_marginalise->d.allele[0]) > 1 || strlen(bcf_record_to_marginalise->d.allele[1]) > 1)
        {
            if (realign(_norm_args, bcf_record_to_marginalise) == ERR_REF_MISMATCH)
            {
                die("vcf record did not match the reference");
            }
        }
        mnp_split(bcf_record_to_marginalise,_hdr,atomised_variants);
    }
    else
    {
        //this block deals with multi-allelics. for each allele a pseudo-unary version is created. with depth and likelihoods marginalsation into a symbolic X allele
        const int num_allele = 3;
        char **new_alleles = new char *[num_allele];
        char *symbolic_allele = "X";
        new_alleles[2] = symbolic_allele;
        int32_t *new_ad = new int32_t[num_allele];
        int32_t *new_gt = new int32_t[2];
        int32_t *new_pl = new int32_t[num_allele * (1 + num_allele) / 2];
        new_ad[0] = old_ad[0];
        for (int i = 1; i < bcf_record_to_marginalise->n_allele; i++)
        {
            new_ad[i]=old_ad[i];
            new_ad[2]=0;
            bcf1_t *new_record = bcf_init1();
            new_record->rid = bcf_record_to_marginalise->rid;
            new_record->pos = bcf_record_to_marginalise->pos;
            bcf_update_format_int32(_hdr, new_record, "GQ", old_gq, num_gq);
            bcf_update_format_int32(_hdr, new_record, "DP", old_dp, num_dp);
            bcf_update_format_int32(_hdr, new_record, "DPF", old_dpf, num_dpf);

            new_alleles[0] = bcf_record_to_marginalise->d.allele[0];
            new_alleles[i] = bcf_record_to_marginalise->d.allele[i];
            bcf_update_alleles(_hdr, new_record, new_alleles);
            for (int j = 0; j < bcf_record_to_marginalise->n_allele; j++)
            {
                if (i != j)//j is an alternate allele that is not i
                {
                    new_ad[2]+=old_ad[j];
                }
            }
            bcf_update_format_int32(_hdr,new_record,"PL",new_pl,num_pl);
            bcf_update_format_int32(_hdr,new_record,"AD",new_ad,num_ad);
            bcf_update_genotypes(_hdr, new_record, new_gt, num_gt);
            if (realign(_norm_args, new_record) == ERR_REF_MISMATCH)
            {
                die("vcf record did not match the reference");
            }
        }
        delete new_alleles;
        delete new_ad;
        delete new_gt;
        delete new_pl;
    }

    delete old_ad;
    delete old_gt;
    delete old_gq;
    delete old_dpf;
    delete old_pl;
    delete old_dp;

    return(atomised_variants);
}


//1. split multi-allelics
//2. normalise (left-align + trim)
//3. decompose MNPs into SNPs
vector<bcf1_t *> Normaliser::atomise(bcf1_t *bcf_record_to_canonicalise)
{
    bcf_unpack(bcf_record_to_canonicalise, BCF_UN_ALL);
    assert(bcf_record_to_canonicalise->n_allele > 1);
    vector<bcf1_t *> atomised_variants;
    bcf1_t **split_records = &bcf_record_to_canonicalise;
    int num_split_records = 1;
    if (bcf_record_to_canonicalise->n_allele > 2)
    {//split multi-allelics (using vcfnorm.c from bcftools1.3
        split_multiallelic_to_biallelics(_norm_args, bcf_record_to_canonicalise);
        split_records = _norm_args->tmp_lines;
        num_split_records = _norm_args->ntmp_lines;
    }

    for (int rec_index = 0; rec_index < num_split_records; rec_index++)
    {
        bcf1_t *rec = split_records[rec_index];

//	std::cerr << "atomise: "<<rec->pos+1<<":"<<rec->d.allele[0]<<":"<<rec->d.allele[1]<<std::endl;
        if (strlen(rec->d.allele[0]) == 1 && strlen(rec->d.allele[1]) == 1)//is a snp. do nothing
        {
            atomised_variants.push_back(bcf_dup(rec));
        }
        else
        {
            if (realign(_norm_args, rec) == ERR_REF_MISMATCH)
            {
                die("vcf record did not match the reference");
            }
            char *ref = rec->d.allele[0];
            char *alt = rec->d.allele[1];
            int ref_len = strlen(ref);
            int alt_len = strlen(alt);

            if (ref_len > 1 && ref_len == alt_len) //is MNP
            {
                char alleles[4] = "X,X";
                for (int i = 0; i < ref_len; i++)
                {
                    if (ref[i] != alt[i])
                    {//new SNP
                        bcf1_t *new_var = bcf_dup(rec);
                        bcf_unpack(new_var, BCF_UN_ALL);
                        alleles[0] = ref[i];
                        alleles[2] = alt[i];
                        new_var->pos += i;
                        bcf_update_alleles_str(_hdr, new_var, alleles);
                        atomised_variants.push_back(new_var);
//			std::cerr << "new_var: "<<new_var->pos+1<<":"<<new_var->d.allele[0]<<":"<<new_var->d.allele[1]<<std::endl;
                    }
                }
            }
            else //variant already is atomic
            {
                bcf1_t *new_var = bcf_dup(rec);
                atomised_variants.push_back(new_var);
            }
        }
    }
    return (atomised_variants);
}
