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
