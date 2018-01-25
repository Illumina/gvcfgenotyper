#include <htslib/vcf.h>
#include "Normaliser.hh"

//#define DEBUG

Normaliser::Normaliser(const string &ref_fname, bool ignore_non_matching_ref) {
    _lg = spdlog::get("gg_logger");
    _norm_args = init_vcfnorm(nullptr, (char *) ref_fname.c_str());
    _symbolic_allele[0] = 'X';
    _symbolic_allele[1] = '\0';
    _ignore_non_matching_ref = ignore_non_matching_ref;
}

// delete/free of symbolic alleles/_hdr
Normaliser::~Normaliser() {
    destroy_data(_norm_args);
    free(_norm_args);
}

int mnp_decompose(bcf1_t *record_to_split, bcf_hdr_t *header, vector<bcf1_t *> &output) {
    int num_allele = record_to_split->n_allele;
    char **alleles = record_to_split->d.allele;
    size_t ref_len = strlen(alleles[0]);
    bool is_mnp = ref_len > 1;
    for (int i = 0; i < num_allele; i++) {
        is_mnp = is_mnp && ref_len == strlen(alleles[i]);
    }

    if (!is_mnp) {
        output.push_back(bcf_dup(record_to_split));
        return (1);
    } else {
        char **new_alleles = (char **) malloc(sizeof(char *) * num_allele);
        for (int i = 0; i < num_allele; i++) {
            new_alleles[i] = (char *) malloc(sizeof(char) * 2);
            new_alleles[i][1] = '\0';
        }
        Genotype old_genotype(header, record_to_split);

        int num_new_snps = 0;
        for (size_t i = 0; i < ref_len; i++) {
            //figures out the set of SNPs that are present at position i in the MNP (might be 0)
            int num_new_allele = 1;
            new_alleles[0][0] = alleles[0][i];
            for (int j = 1; j < num_allele; j++) {
                if (alleles[0][i] != alleles[j][i]) {
                    bool has_been_seen = false;
                    for (int k = 0; k < j; k++) {
                        has_been_seen |= alleles[j][i] == alleles[k][i];
                    }
                    if (!has_been_seen) {
                        num_new_allele++;
                        new_alleles[num_new_allele - 1][0] = alleles[j][i];
                    }
                }
            }

            //if there are new alternate alleles present, this manages the FORMAT fields
            if (num_new_allele > 1) {
                bcf1_t *new_var = bcf_dup(record_to_split);
                bcf_unpack(new_var, BCF_UN_ALL);
                new_var->pos += i;
                bcf_update_alleles(header, new_var, (const char **) new_alleles, num_new_allele);
                //the number of alleles changed so we have to reformat the FORMAT fields
                if (num_new_allele != num_allele) {
                    //creates a mapping from old alleles -> new alleles
                    Genotype new_genotype(old_genotype._ploidy, num_new_allele);
                    new_genotype.SetDepthToZero();
                    vector<int> allele_remap(num_allele);
                    for (int new_allele = 0; new_allele < num_new_allele; new_allele++) {
                        for (int old_allele = 0; old_allele < num_allele; old_allele++) {
                            if (alleles[old_allele][i] == new_alleles[new_allele][0]) {
                                allele_remap[old_allele] = new_allele;
                            }
                        }
                    }
                    //remaps old genotypes to new genotypes
                    for (int genotype_index = 0; genotype_index < old_genotype._ploidy; genotype_index++) {
                        if (bcf_gt_is_phased(old_genotype._gt[genotype_index])) {
                            new_genotype._gt[genotype_index] = bcf_gt_phased(
                                    allele_remap[bcf_gt_allele(old_genotype._gt[genotype_index])]);
                        } else {
                            new_genotype._gt[genotype_index] = bcf_gt_unphased(
                                    allele_remap[bcf_gt_allele(old_genotype._gt[genotype_index])]);
                        }
                    }
                    //marginalises FORMAT/AD
                    for (int i = 0; i < num_allele; i++) {
                        new_genotype._ad[allele_remap[i]] += old_genotype._ad[i];
                    }
                    //marginalises FORMAT/PL
                    for (int i = 0; i < num_allele; i++) {
                        for (int j = i; j < num_allele; j++) {
                            new_genotype._gl[ggutils::get_gl_index(allele_remap[i],
                                                                   allele_remap[j])] += old_genotype._gl[ggutils::get_gl_index(
                                    i, j)];
                        }
                    }
                    new_genotype.SetDp(old_genotype.dp());
                    new_genotype.SetDpf(old_genotype.dpf());
                    new_genotype.SetGq(old_genotype.gq());
                    new_genotype.SetGqx(old_genotype.gqx());
                    new_genotype.UpdateBcfRecord(header, new_var);
                }
                output.push_back(new_var);
                num_new_snps++;
            }
        }
        for (int i = 0; i < num_allele; i++) {
            free(new_alleles[i]);
        }
        free(new_alleles);
        return (num_new_snps);
    }
}

//Performs left-alignment and trimming using code from bcftools' vcfnorm.c
bool Normaliser::Realign(bcf1_t *record, bcf_hdr_t *header) {
    if (realign(_norm_args, record, header) != ERR_OK) {
        if (_ignore_non_matching_ref) {
            _lg->warn("WARNING: VCF record did not match the reference at sample {}", header->samples[0]);
            return (false);
        } else {
            ggutils::die("VCF record did not match the reference at sample " + (string) header->samples[0]);
        }
    }
    return (true);
}

void Normaliser::MultiSplit(bcf1_t *bcf_record_to_split, vector<bcf1_t *> &split_variants, bcf_hdr_t *hdr) {
    assert(bcf_record_to_split->n_allele > 2);
    bcf_unpack(bcf_record_to_split, BCF_UN_ALL);
    Genotype src(hdr, bcf_record_to_split);
    Genotype dst(src.ploidy(), src.num_allele());

    std::vector<std::pair<int, int> > new_positions; //stores the position + rank of each variant post-normalisation
    char **new_alleles = (char **) malloc(sizeof(char *) * bcf_record_to_split->n_allele);
    for (int i = 1; i < bcf_record_to_split->n_allele; i++) {
        bcf1_t *tmp_record = bcf_dup(bcf_record_to_split);
        bcf_unpack(tmp_record, BCF_UN_ALL);
        new_alleles[0] = bcf_record_to_split->d.allele[0];
        new_alleles[1] = bcf_record_to_split->d.allele[i];
        bcf_update_alleles(hdr, tmp_record, (const char **) new_alleles, 2);
        if (Realign(tmp_record, hdr))
            new_positions.push_back(pair<int, int>(tmp_record->pos, ggutils::get_variant_rank(tmp_record)));
        bcf_destroy(tmp_record);
    }

    std::set<std::pair<int, int> > unique_positions(new_positions.begin(), new_positions.end());

    for (auto pos = unique_positions.begin(); pos != unique_positions.end(); pos++) {
        int counter = 1;
        vector<int> alleles_at_this_position;
        new_alleles[0] = bcf_record_to_split->d.allele[0];
        for (size_t i = 0; i < new_positions.size(); i++) {
            if (new_positions[i] == *pos) {
                alleles_at_this_position.push_back(1 + (int) i);
                new_alleles[counter++] = bcf_record_to_split->d.allele[i + 1];
            }
        }
        assert(alleles_at_this_position.size() > 0);

        bcf1_t *tmp_record = bcf_dup(bcf_record_to_split);
        bcf_unpack(tmp_record, BCF_UN_ALL);
        bcf_update_alleles(hdr, tmp_record, (const char **) new_alleles, 1 + (int) alleles_at_this_position.size());

        if (unique_positions.size() > 1) {
            src.CollapseAllelesIntoRef(alleles_at_this_position, dst);
            dst.UpdateBcfRecord(hdr, tmp_record);
        } else {
            src.UpdateBcfRecord(hdr, tmp_record);
        }
        for (int i = 1; i < tmp_record->n_allele; i++) {
            bcf1_t *out_record = bcf_dup(tmp_record);
            bcf_unpack(out_record, BCF_UN_ALL);
            ggutils::bcf1_allele_swap(hdr, out_record, i, 1);
            if (Realign(out_record, hdr))split_variants.push_back(out_record);
        }
        bcf_destroy1(tmp_record);
    }
    free(new_alleles);
}

void Normaliser::Unarise(bcf1_t *bcf_record_to_marginalise, vector<bcf1_t *> &atomised_variants, bcf_hdr_t *hdr) {
#ifdef DEBUG
    ggutils::print_variant(hdr,bcf_record_to_marginalise);
#endif
    //bi-allelic snp. Nothing to do, just copy the variant into the buffer.
    if (ggutils::is_snp(bcf_record_to_marginalise) && bcf_record_to_marginalise->n_allele == 2) {
        atomised_variants.push_back(bcf_dup(bcf_record_to_marginalise));
        return;
    }

    //FIXME: We would like to get rid of this special-case MNP decomposition and replace it with a more general decomposition step.
    //FIXME: For now this at least allows us to behave well for SNPs that are hidden in MNPS.
    vector<bcf1_t *> decomposed_variants;
    mnp_decompose(bcf_record_to_marginalise, hdr, decomposed_variants);

    for (auto it = decomposed_variants.begin(); it != decomposed_variants.end(); ++it) {
        bcf1_t *decomposed_record = *it;
        if (decomposed_record->n_allele == 2)//bi-allelic. no further decomposition needed.
        {
            if (Realign(decomposed_record, hdr))
                atomised_variants.push_back(decomposed_record);
        } else {
            MultiSplit(*it, atomised_variants, hdr);
            bcf_destroy(*it);
        }
    }
}

bcf1_t *CollapseRecords(bcf_hdr_t *sample_header,
                        pair<std::deque<bcf1_t *>::iterator, std::deque<bcf1_t *>::iterator> &sample_variants) {
    if ((sample_variants.second - sample_variants.first) == 0)
        return nullptr;

    bcf1_t *ret = bcf_dup(*sample_variants.first);
    bcf_unpack(ret, BCF_UN_ALL);
    if (sample_variants.first == sample_variants.second)
        return (ret);

    int num_old_alleles = ret->n_allele;
    int num_new_alleles = 0;
    for (auto it = (sample_variants.first + 1); it != sample_variants.second; it++)
        for (int i = 1; i < (*it)->n_allele; i++)
            num_new_alleles += ggutils::add_allele(sample_header, ret, *it, i) > num_old_alleles;

    if (num_new_alleles == 0)
        return (ret);

    //This is where things get messy. We have a collision between overlapping alleles that
    //do not have respsective INFO/FORMAT for one another.
    int ploidy = 0;
    size_t num_allele = (sample_variants.second - sample_variants.first) + 1;
    for (auto it = sample_variants.first; it != sample_variants.second; it++)
        ploidy = max(ggutils::get_ploidy(sample_header, *it), ploidy);
    assert(ploidy == 1 || ploidy == 2);
    Genotype output(ploidy, num_allele);
    std::cerr << "num_pl"<<output.num_pl() << std::endl;//debug
    std::vector<bool> found_allele(sample_variants.second - sample_variants.first + 1, false);
    std::vector<std::vector<int> > pls;
    for (auto it = sample_variants.first; it != sample_variants.second; it++) {
        int allele_index = ggutils::find_allele(ret, *it, 1);
        if (allele_index == -1 || !found_allele[allele_index]) {
            Genotype g(sample_header, *it);
            pls.emplace_back(output.num_pl(), bcf_int32_missing);
            for (int i = 0; i < (*it)->n_allele; i++) {
                allele_index = i == 0 ? 0 : ggutils::find_allele(ret, *it, i);
                found_allele[allele_index] = true;
                if (output.ad(allele_index) == bcf_int32_missing) {
                    output.SetAd(g.ad(i), allele_index);
                    output.SetAdf(g.adf(i), allele_index);
                    output.SetAdr(g.adr(i), allele_index);
                }
                for (int j = 0; j < (*it)->n_allele; j++) {
                    int allele_index2 = j == 0 ? 0 : ggutils::find_allele(ret, *it, j);
                    pls.back()[ggutils::get_gl_index(allele_index, allele_index2)] = g.pl(i, j);
                }
            }
        }
    }
    std::vector<int> pl;
    ggutils::collapse_gls(ploidy, num_allele, pls, pl);
    output.SetPl(pl);
    output.SetDepthFromAd();
    output.CallGenotype();
    output.UpdateBcfRecord(sample_header, ret);
    return (ret);
}
