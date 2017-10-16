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
    int num_allele = record_to_split->n_allele;
    char **alleles = record_to_split->d.allele;
    int ref_len = strlen(alleles[0]);
    bool is_mnp = ref_len>1;
    char **new_alleles = new char *[num_allele];
    for(int i=0;i<num_allele;i++)
    {
        new_alleles[i] = new char[2];
        new_alleles[i][1] = '\0';
        is_mnp = is_mnp && ref_len==strlen(alleles[i]);
    }

    int num_new_snps=0;
    if(is_mnp)
    {
        int *new_gt=NULL,nval=0;
        int ploidy = bcf_get_genotypes(header,record_to_split,&new_gt,&nval);
        int num_gl = ploidy==1 ? ploidy : num_allele * (1 + num_allele) / 2;
        auto old_ad = new int32_t[num_allele];
        auto new_ad = new int32_t[num_allele];
        auto old_pl = new int32_t[num_gl];
        auto new_pl = new int32_t[num_gl];
        assign(num_allele,0,old_ad);
        assign(num_allele,0,new_ad);
        assign(num_gl,0,old_pl);
        assign(num_gl,0,new_pl);
        assert(bcf_get_format_int32(header,record_to_split,"AD",&old_ad, &num_allele) == num_allele );
        assert(bcf_get_format_int32(header,record_to_split,"PL",&old_pl, &num_gl) == num_gl );
        for (int i = 0; i < ref_len; i++)
        {
            int num_new_alts = 0;
            new_alleles[0][0] = alleles[0][i];
            for(int j=1;j<num_allele;j++)
            {
                if(alleles[0][i]!=alleles[j][i])
                {
                    bool has_been_seen = false;
                    for(int k=0;k<j;k++)
                    {
                        has_been_seen = alleles[j][i]==alleles[k][i];
                    }
                    if(!has_been_seen)
                    {
                        num_new_alts++;
                        new_alleles[num_new_alts][0] = alleles[j][i];
                    }
                }
            }

            if(num_new_alts>0)
            {
                int num_new_gl = ploidy==1 ? ploidy : (num_new_alts+1) * (2 + num_new_alts) / 2;
                bcf1_t *new_var = bcf_dup(record_to_split);
                bcf_unpack(new_var, BCF_UN_ALL);
                new_var->pos += i;
                bcf_update_alleles(header, new_var, (const char **)new_alleles,num_new_alts+1);
                ploidy = bcf_get_genotypes(header,record_to_split,&new_gt,&nval);

                for(int new_allele=0;new_allele<num_new_alts;new_allele++)
                {
                    for(int old_allele=0;old_allele<num_allele;old_allele++)
                    {
                        if(alleles[old_allele][i]==new_alleles[new_allele][0])
                        {
                            for(int genotype_index=0;genotype_index<ploidy;genotype_index++)
                            {
                                if(bcf_gt_allele(new_gt[genotype_index])==old_allele)
                                {
                                    if(bcf_gt_is_phased(new_gt[genotype_index]))
                                    {
                                        new_gt[genotype_index] = bcf_gt_phased(new_allele);
                                    }
                                    else
                                    {
                                        new_gt[genotype_index] = bcf_gt_unphased(new_allele);
                                    }
                                }
                            }
                        }
                    }
                }
                bcf_update_genotypes(header,new_var,new_gt,ploidy);
                bcf_update_format_int32(header,new_var,"AD",&(new_ad),num_new_alts);
                bcf_update_format_int32(header,new_var,"PL",&(new_pl),num_new_gl);
                output.push_back(new_var);
                num_new_snps++;
            }
        }
        for(int i=0;i<num_allele;i++)
        {
            delete new_alleles[i];
        }
        delete[] new_alleles;
        return(num_new_snps);
    }
    else
    {
        output.push_back(bcf_dup(record_to_split));
        return(1);
    }
}

static int factorial(int x)
{
    int ret = 1;
    for(int i=2;i<=x;i++)
    {
        ret *= i;
    }
    return(ret);
}

static int choose(int n,int k)
{
    return(factorial(n)/(factorial(n-k)*factorial(k)));
}

static int get_genotype_index(int g0,int g1)
{
    return(choose(g0 + 1 - 1,1) + choose(g1 + 2 - 1,2));
}
//
//vector<bcf1_t *> Normaliser::unarise(bcf1_t *bcf_record_to_marginalise)
//{
//    int old_num_allele = bcf_record_to_marginalise->n_allele;
//    bcf_unpack(bcf_record_to_marginalise, BCF_UN_ALL);
//    assert(old_num_allele > 1);
//    vector<bcf1_t *> atomised_variants;
//
//    //this chunk of codes reads the relevant FORMAT fields (PL,GQ,DP,DPF,AD) that will be propagated into the new record.
//    int32_t *old_ad,*old_gt,*old_gq,*old_dpf,*old_pl,*old_dp;
//    int ploidy = bcf_get_genotypes(_hdr,bcf_record_to_marginalise,&old_gt,&num_gt);
//    assert(ploidy>=0 && ploidy<=2);
//
//    int num_gl = ploidy==1 ? ploidy : old_num_allele * (1 + old_num_allele) / 2;
//    int num_ad=0,num_gt=0,num_gq=0,num_dpf=0,num_pl=0,num_dp=0;
//    assert(bcf_get_format_int32(_hdr,bcf_record_to_marginalise,"PL",&old_pl,&num_pl)==num_gl);
//    assert(bcf_get_format_int32(_hdr,bcf_record_to_marginalise,"AD",&old_ad,&num_ad)==old_num_allele);
//    assert(bcf_get_format_int32(_hdr,bcf_record_to_marginalise,"DPF",&old_dpf,&num_dpf)==1);
//    assert(bcf_get_format_int32(_hdr,bcf_record_to_marginalise,"DP",&old_dp,&num_dp)==1);
//    if(bcf_get_format_int32(_hdr,bcf_record_to_marginalise,"GQ",&old_gq,&num_gq)==-2)//catches cases where GQ is set as a float (should be int according to vcf4.3)
//    {
//        float *tmp_gq;
//        assert(bcf_get_format_float(_hdr,bcf_record_to_marginalise,"GQ",&tmp_gq,&num_gq)==1);
//        old_gq = new int32_t[1];
//        old_gq[0] = (int32_t)tmp_gq[0];
//        delete tmp_gq;
//    }
//
//    //this block deals with bi-allelic variants. no marginalisation across other alternate alleles is required
//    if(old_num_allele == 2)
//    {
//        bcf1_t *new_record = bcf_init1();
//        new_record->rid = bcf_record_to_marginalise->rid;
//        new_record->pos = bcf_record_to_marginalise->pos;
//        bcf_update_alleles(_hdr,new_record,bcf_record_to_marginalise->d.allele,old_num_allele);
//        bcf_update_genotypes(_hdr,new_record,old_gt,num_gt);
//        bcf_update_format_int32(_hdr,new_record,"GQ",old_gq,num_gq);
//        bcf_update_format_int32(_hdr,new_record,"AD",old_ad,num_ad);
//        bcf_update_format_int32(_hdr,new_record,"DP",old_dp,num_dp);
//        bcf_update_format_int32(_hdr,new_record,"DPF",old_dpf,num_dpf);
//        bcf_update_format_int32(_hdr,new_record,"PL",old_pl,num_pl);
//
//        if(strlen(bcf_record_to_marginalise->d.allele[0]) > 1 || strlen(bcf_record_to_marginalise->d.allele[1]) > 1)
//        {
//            if (realign(_norm_args, bcf_record_to_marginalise) == ERR_REF_MISMATCH)
//            {
//                die("vcf record did not match the reference");
//            }
//        }
//        mnp_split(bcf_record_to_marginalise,_hdr,atomised_variants);
//    }
//    else
//    {
//        //this block deals with multi-allelics. for each allele a pseudo-unary version is created. with depth and likelihoods marginalsation into a symbolic X allele
//        const int num_allele = 3;
//        auto **new_alleles = new char *[num_allele];
//        char *symbolic_allele = "X";
//        new_alleles[2] = symbolic_allele;
//        auto *new_ad = new int32_t[num_allele];
//        auto *new_gt = new int32_t[2];
//        auto *new_pl = new int32_t[num_gl];
//        vector<float> old_gl(num_gl,0.);
//        vector<float> new_gl(num_gl / 2,0.);
//        for(int i=0;i<num_gl;i++)
//        {
//            old_gl[i] = pow((float)10.,(float)-old_pl[i]/10.);
//        }
//        int reference_allele=0;
//        int primary_allele=1;
//        int symbolic_allele=2;
//        new_ad[reference_allele] = old_ad[reference_allele];
//        for (int i = 1; i<old_num_allele ; i++)
//        {
//            new_ad[primary_allele]=old_ad[i];
//            new_ad[symbolic_allele]=0;
//            bcf1_t *new_record = bcf_init1();
//            new_record->rid = bcf_record_to_marginalise->rid;
//            new_record->pos = bcf_record_to_marginalise->pos;
//            bcf_update_format_int32(_hdr, new_record, "GQ", old_gq, num_gq);
//            bcf_update_format_int32(_hdr, new_record, "DP", old_dp, num_dp);
//            bcf_update_format_int32(_hdr, new_record, "DPF", old_dpf, num_dpf);
//
//            new_alleles[0] = bcf_record_to_marginalise->d.allele[0];
//            new_alleles[i] = bcf_record_to_marginalise->d.allele[i];
//            bcf_update_alleles(_hdr, new_record, new_alleles);
//            for(int j=0;j<ploidy;j++)
//            {
//                if(bcf_gt_allele(old_gt[j])==reference_allele)
//                {
//                    new_gt[j] = bcf_gt_unphased(reference_allele);
//                }
//                else if(bcf_gt_allele(old_gt[j])==i)
//                {
//                    new_gt[j] = bcf_gt_unphased(primary_allele);
//                }
//                else
//                {
//                    new_gt[j] = bcf_gt_unphased(symbolic_allele);
//                }
//            }
//
//            //marginalises FORMAT/AD
//            for (int j = 1; j < old_num_allele; j++)
//            {
//                if (i != j)//j is an alternate allele that is not i
//                {
//                    new_ad[symbolic_allele]+=old_ad[j];
//                }
//            }
//            bcf_update_format_int32(_hdr,new_record,"AD",new_ad,num_ad);
//
//            //marginalises FORMAT/PL
////            new_gl.assign(num_gl,0.);
////            int count = 0;
////            for(int j=0;j<old_num_allele;j++)
////            {
////                if(ploidy==1)
////                {
////                    if(j==reference_allele)
////                    {
////                        new_gl[reference_allele] = old_gl[j];
////                    }
////                    else if(i==j)
////                    {
////                        new_gl[primary_allele] = old_gl[j];
////                    }
////                    else
////                    {
////                        new_gl[symbolic_allele] += old_gl[j];
////                    }
////                }
////                if(ploidy==2)
////                {
////                    for(int k=0;k<ploidy;k++)
////                    {
////                        if(j==reference_allele)
////                        {
////                            if(k==reference_allele)
////                            {
////                                new_gl[get_genotype_index(reference_allele,reference_allele)] +=  old_gl[count];
////                            }
////                            else if(k==i)
////                            {
////                                new_gl[get_genotype_index(reference_allele,primary_allele)] +=  old_gl[count];
////                            }
////                            else
////                            {
////                                new_gl[get_genotype_index(reference_allele,symbolic_allele)] +=  old_gl[count];
////                            }
////                        }
////                        else if()
////                        count++;
////                    }
////                }
////            }
//            new_gl.assign(num_gl,1.);
//            for(int j=0;j<num_gl;j++)
//            {
//                new_pl[j]=-(int32_t)log10(10*new_gl[j]);
//            }
//            bcf_update_format_int32(_hdr,new_record,"PL",new_pl,num_pl);
//
//            bcf_update_genotypes(_hdr, new_record, new_gt, num_gt);
//            if (realign(_norm_args, new_record) == ERR_REF_MISMATCH)
//            {
//                die("vcf record did not match the reference");
//            }
//        }
//        delete new_alleles;
//        delete new_ad;
//        delete new_gt;
//        delete new_pl;
//    }
//
//    delete old_ad;
//    delete old_gt;
//    delete old_gq;
//    delete old_dpf;
//    delete old_pl;
//    delete old_dp;
//
//    return(atomised_variants);
//}
//

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
