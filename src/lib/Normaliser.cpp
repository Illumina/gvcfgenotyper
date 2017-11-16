#include "Normaliser.hh"

// remove all INFO fields
/*static void remove_info(bcf1_t *line)
{
    if (!(line->unpacked & BCF_UN_INFO))
    { 
        bcf_unpack(line, BCF_UN_INFO); 
    }
    
    for (int i = 0; i < line->n_info; i++)
    {
        bcf_info_t *inf = &line->d.info[i];
        if (inf->vptr_free)
        {
            free(inf->vptr - inf->vptr_off);
            inf->vptr_free = 0;
        }
        line->d.shared_dirty |= BCF1_DIRTY_INF;
        inf->vptr = nullptr;
    }
    line->n_info = 0;
}*/

Normaliser::Normaliser(const string &ref_fname, bcf_hdr_t *hdr)
{
    _hdr = hdr;
    _norm_args = init_vcfnorm(_hdr, (char *)ref_fname.c_str());
    _symbolic_allele[0]='X';
    _symbolic_allele[1]='\0';
}

// delete/free of symbolic alleles/_hdr
Normaliser::~Normaliser()
{
    destroy_data(_norm_args);
    free(_norm_args);
}

int mnp_split(bcf1_t *record_to_split, bcf_hdr_t *header, vector<bcf1_t *> & output)
{
    int num_allele = record_to_split->n_allele;
    char **alleles = record_to_split->d.allele;
    size_t ref_len = strlen(alleles[0]);
    bool is_mnp = ref_len > 1;
    for (int i = 0; i < num_allele; i++)
    {
        is_mnp = is_mnp && ref_len == strlen(alleles[i]);
    }

    if (!is_mnp)
    {
        output.push_back(bcf_dup(record_to_split));
        return (1);
    }
    else
    {
        char **new_alleles = (char **)malloc(sizeof(char *)*num_allele);
        for (int i = 0; i < num_allele; i++)
        {
            new_alleles[i] = (char *)malloc(sizeof(char)*2);
            new_alleles[i][1] = '\0';
        }
        Genotype old_genotype(header, record_to_split);

        int num_new_snps = 0;
        for (size_t i = 0; i < ref_len; i++)
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
                //the number of alleles changed so we have to reformat the FORMAT fields
                if (num_new_allele != num_allele)
                {
                    //creates a mapping from old alleles -> new alleles
                    Genotype new_genotype(old_genotype._ploidy, num_new_allele);
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
                    new_genotype._dp[0] = old_genotype._dp[0];
                    new_genotype._dpf[0] = old_genotype._dpf[0];
                    new_genotype._gq[0] = old_genotype._gq[0];
                    new_genotype.update_bcf1_t(header, new_var);
                }
                output.push_back(new_var);
                num_new_snps++;
            }
        }
        for (int i = 0; i < num_allele; i++)
        {
            free(new_alleles[i]);
        }
        free(new_alleles);
        return (num_new_snps);
    }
}


void Normaliser::unarise(bcf1_t *bcf_record_to_marginalise, vector<bcf1_t*>& atomised_variants )
{
    //bi-allelic snp. Nothing to do, just copy the variant into the buffer.
    if(is_snp(bcf_record_to_marginalise) && bcf_record_to_marginalise->n_allele==2)
    {
        atomised_variants.push_back(bcf_dup(bcf_record_to_marginalise));
        return;
    }

    const int reference_allele = 0;
    const int primary_allele = 1;
    const int symbolic_allele = 2;
    const int num_new_allele = 3;

    auto **new_alleles = new char *[num_new_allele];

    //FIXME: we would like to get rid of this special-case MNP decomposition and replace it with a more general decomposition step.
    //FIXME: for now this at least allows us to behave well for SNPs that are hidden in MNPS
    vector<bcf1_t *> decomposed_variants;
    mnp_split(bcf_record_to_marginalise, _hdr, decomposed_variants);

    for (auto it = decomposed_variants.begin(); it != decomposed_variants.end(); ++it)
    {
        bcf1_t *decomposed_record = *it;
        if(decomposed_record->n_allele==2)//bi-alleic. no further decomposition needed.
        {
            if (realign(_norm_args, decomposed_record) != ERR_OK)
            {
                die("vcf record did not match the reference");
            }
            atomised_variants.push_back(decomposed_record);
        }
        else
        {
            Genotype old_genotype(_hdr, decomposed_record);
            for (int i = 1; i < old_genotype._num_allele; i++)
            {
                bcf1_t *new_record = bcf_dup(decomposed_record);
                bcf_unpack(new_record,BCF_UN_ALL);
                new_alleles[reference_allele] = decomposed_record->d.allele[reference_allele];
                new_alleles[primary_allele] = decomposed_record->d.allele[i];
                bcf_update_alleles(_hdr, new_record, (const char **) new_alleles, num_new_allele-1);
                if (realign(_norm_args, new_record) != ERR_OK)
                {
                    die("vcf record did not match the reference");
                }
                //now add the symbolic allele
                new_alleles[reference_allele] =  new_record->d.allele[0];
                new_alleles[primary_allele] =  new_record->d.allele[1];
                new_alleles[symbolic_allele] =  _symbolic_allele;
                bcf_update_alleles(_hdr, new_record, (const char **) new_alleles, num_new_allele);

                Genotype new_genotype = old_genotype.marginalise(i);
                new_genotype.update_bcf1_t(_hdr, new_record);
                atomised_variants.push_back(new_record);
            }
            bcf_destroy(decomposed_record);
        }
    }
    delete[] new_alleles;
}
