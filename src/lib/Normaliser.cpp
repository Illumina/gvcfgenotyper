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

int mnp_split(bcf1_t *record_to_split, bcf_hdr_t *header, vector<bcf1_t *> &output)
{
    int num_allele = record_to_split->n_allele;
    char **alleles = record_to_split->d.allele;
    int ref_len = strlen(alleles[0]);
    bool is_mnp = ref_len > 1;
    for (int i = 0; i < num_allele; i++)
    {
        is_mnp = is_mnp && ref_len == strlen(alleles[i]);
    }

    if (is_mnp)
    {
        char **new_alleles = new char *[num_allele];
        for (int i = 0; i < num_allele; i++)
        {
            new_alleles[i] = new char[2];
            new_alleles[i][1] = '\0';
        }
        Genotype old_genotype(header, record_to_split);
        int ploidy = old_genotype._ploidy;
        int num_new_snps = 0;
        for (int i = 0; i < ref_len; i++)
        {
            //figures out the set of SNPs that are present at position i in the MNP (might be 0)
            int num_new_allele = 1;
            new_alleles[0][0] = alleles[0][i];
            for (int j = 1; j < num_allele; j++)
            {
                if (alleles[0][i] != alleles[j][i])
                {
                    bool has_been_seen = false;
                    for (int k = 0; k < j; k++)
                    {
                        has_been_seen |= alleles[j][i] == alleles[k][i];
                    }
                    if (!has_been_seen)
                    {
                        num_new_allele++;
                        new_alleles[num_new_allele - 1][0] = alleles[j][i];
                    }
                }
            }

            //if there are new alternate alleles present, this manages the FORMAT fields
            if (num_new_allele > 1)
            {
                bcf1_t *new_var = bcf_dup(record_to_split);
                bcf_unpack(new_var, BCF_UN_ALL);
                new_var->pos += i;
                bcf_update_alleles(header, new_var, (const char **) new_alleles, num_new_allele);
                if (num_new_allele !=
                    num_allele)//the number of alleles changed so we have to reformat the FORMAT fields
                {
                    Genotype new_genotype(old_genotype._ploidy, num_new_allele);

                    //creates a mapping from old alleles -> new alleles
                    vector<int> allele_remap(num_allele);
                    for (int new_allele = 0; new_allele < num_new_allele; new_allele++)
                    {
                        for (int old_allele = 0; old_allele < num_allele; old_allele++)
                        {
                            if (alleles[old_allele][i] == new_alleles[new_allele][0])
                            {
                                allele_remap[old_allele] = new_allele;
                            }
                        }
                    }
                    //remaps old genotypes to new genotypes
                    for (int genotype_index = 0; genotype_index < old_genotype._ploidy; genotype_index++)
                    {
                        if (bcf_gt_is_phased(old_genotype._gt[genotype_index]))
                        {
                            new_genotype._gt[genotype_index] = bcf_gt_phased(
                                    allele_remap[bcf_gt_allele(old_genotype._gt[genotype_index])]);
                        }
                        else
                        {
                            new_genotype._gt[genotype_index] = bcf_gt_unphased(
                                    allele_remap[bcf_gt_allele(old_genotype._gt[genotype_index])]);
                        }
                    }
                    //marginalises FORMAT/AD
                    for (int i = 0; i < num_allele; i++)
                    {
                        new_genotype._ad[allele_remap[i]] += old_genotype._ad[i];
                    }
                    //marginalises FORMAT/PL
                    for (int i = 0; i < num_allele; i++)
                    {
                        for (int j = i; j < num_allele; j++)
                        {
                            new_genotype._gl[get_gl_index(allele_remap[i],
                                                          allele_remap[j])] += old_genotype._gl[get_gl_index(i, j)];
                        }
                    }
                    new_genotype._dpf[0] = old_genotype._dpf[0];
                    new_genotype._gq[0] = old_genotype._gq[0];
                    new_genotype.setDepthFromAD();
                    new_genotype.update_bcf1_t(header, new_var);
                }
                output.push_back(new_var);
                num_new_snps++;
            }
        }
        for (int i = 0; i < num_allele; i++)
        {
            delete new_alleles[i];
        }
        delete[] new_alleles;
        return (num_new_snps);
    }
    else
    {
        output.push_back(bcf_dup(record_to_split));
        return (1);
    }
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
            mnp_split(rec, _hdr, atomised_variants);
        }
    }
    return (atomised_variants);
}

