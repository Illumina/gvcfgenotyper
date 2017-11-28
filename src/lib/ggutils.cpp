#include "ggutils.hh"

namespace ggutils
{
    void right_trim(const char *ref,const char *alt,size_t &reflen,size_t &altlen)
    {
        reflen=strlen(ref);
        altlen=strlen(alt);
        assert(reflen>0);
        assert(altlen>0);
        while(ref[reflen-1]==alt[altlen-1] && reflen>1 && altlen>1)
        {
            reflen--;
            altlen--;
        }
    }

    int *zeros(int n)
    {
        int *ret = (int *) malloc(sizeof(int) * n);
        for (int i = 0; i < n; i++)
        {
            ret[i] = 0;
        }
        return (ret);
    }

    bool fileexists(const string &fname)
    {
        ifstream ifile(fname.c_str());
        if (ifile)
        {
            return (true);
        }
        else
        {
            return (false);
        }
    }

    string percent(int num, int den)
    {
        char buffer[100];
        sprintf(buffer, "%d / %d\t( %.4f%% )\t", num, den, 100. * float(num) / float(den));
        return (buffer);
    }

    int strsplit(const string &input, const char split, vector<string> &out)
    {
        istringstream ss(input);
        out.clear();
        string tmp;
        int count = 0;
        while (std::getline(ss, tmp, split))
        {
            out.push_back(tmp);
            count++;
        }
        return (count);
    }

    string join(const vector<string> &input, const string &delim)
    {
        string ret = input[0];
        for (size_t i = 1; i < input.size(); i++)
        {
            ret += delim;
            ret += input[i];
        }
        return (ret);
    }

    vector<int> match(const vector<string> &x, const vector<string> &y)
    {
        vector<int> ret;
        for (size_t i = 0; i < y.size(); i++)
        {
            size_t j = 0;
            while (x[j] != y[i])
            {
                j++;
                if (j >= x.size())
                {
                    die("sample in y is not in x");
                }
            }
            ret.push_back(j);
        }
        return (ret);
    }

    //returns the ##contig=<ID=chr19,length=59128983> lines from a bcf/vcf.
    int copy_contigs(const bcf_hdr_t *src, bcf_hdr_t *dst)
    {
        vector<string> parseme;
        kstring_t splitme = {0, 0, nullptr};
        bcf_hdr_format(src, 1, &splitme);
        strsplit(splitme.s, '\n', parseme);
        string contigs = "";
        for (size_t i = 0; i < parseme.size(); i++)
        {
            if (parseme[i].find("contig") < parseme[i].size())
            {
                bcf_hdr_append(dst, parseme[i].c_str());
            }
        }
        free(splitme.s);
        return (0);
    }

    int read_text_file(const string &fname, vector<string> &output)
    {
        ifstream ifile(fname.c_str());
        if (!ifile)
        {
            die("problem opening " + ((string) fname));
        }
        string tmp;
        output.clear();
        while (getline(ifile, tmp))
        {
            output.push_back(tmp);
        }

        return (output.size());
    }

    bool is_snp(bcf1_t *record)
    {
        assert(record->n_allele > 1);
        bcf_unpack(record, BCF_UN_ALL);
        return (bcf_get_variant_type(record, 1) & VCF_SNP);
    }

    bool is_hom_ref(const bcf_hdr_t * header, bcf1_t* record)
    {
        int *gt = nullptr, ngt = 0;
        bool is_hom_ref = false;
        if(bcf_get_genotypes(header, record, &gt, &ngt)<=0)
        {
            // A vcf record without GT cannot be hom ref
            return is_hom_ref;
        }
        int ploidy = bcf_get_genotypes(header, record, &gt, &ngt);
        // check for 0/0
        is_hom_ref = (bcf_gt_allele(gt[0]) == 0);
        is_hom_ref &= ploidy==1 || (bcf_gt_allele(gt[1]) == 0);
        free(gt);
        return is_hom_ref;
    }

    int get_end_of_gvcf_block(bcf_hdr_t *header, bcf1_t *record)
    {
        int ret;
        int *ptr = NULL, nval = 0;
        if (bcf_get_info_int32(header, record, "END", &ptr, &nval) == 1)
        {
            ret = *ptr - 1;
            free(ptr);
        }
        else
        {
            ret = record->pos + strlen(record->d.allele[0]) - 1;
        }

        return (ret);
    }

    int get_end_of_variant(bcf1_t *record)
    {
        return (record->pos + strlen(record->d.allele[0]) - 1);
    }

    bool is_deletion(bcf1_t *record)
    {
        bcf_unpack(record, BCF_UN_ALL);
        int l1 = strlen(record->d.allele[0]);
        int l2 = strlen(record->d.allele[1]);
        return (bcf_get_variant_type(record, 1) == VCF_INDEL && l2 < l1);
    }

    bool is_insertion(bcf1_t *record)
    {
        bcf_unpack(record, BCF_UN_ALL);
        int l1 = strlen(record->d.allele[0]);
        int l2 = strlen(record->d.allele[1]);
        return (bcf_get_variant_type(record, 1) == VCF_INDEL && l2 > l1);
    }

    bool is_complex(bcf1_t *record)
    {
        return (!is_snp(record) && !is_deletion(record) && !is_insertion(record));
    }

    int get_variant_rank(bcf1_t *record)
    {
        if (is_complex(record) || is_snp(record))
        {
            return (0);
        }
        if (is_insertion(record))
        {
            return (1);
        }
        if (is_deletion(record))
        {
            return (2);
        }
        die("bad variant");
        return (-1);
    }

    bool bcf1_equal(bcf1_t *a, bcf1_t *b)
    {
        assert(a->n_allele>1);
        assert(b->n_allele>1);
        bcf_unpack(a, BCF_UN_ALL);
        bcf_unpack(b, BCF_UN_ALL);
        if (a == nullptr || b == nullptr)
        {
            die(" (bcf1_equal: tried to compare NULL bcf1_t");
        }
        if (a->rid != b->rid)
        {
            return (false);
        }
        else if (a->pos != b->pos)
        {
            return (false);
        }
        else
        {
            size_t a1,b1,a2,b2;
            right_trim(a->d.allele[0],a->d.allele[1],a1,b1);
            right_trim(b->d.allele[0],b->d.allele[1],a2,b2);
            if(a1!=a2 || b1!=b2) return(false);
            return(strncmp(a->d.allele[0],b->d.allele[0],a1)==0 && strncmp(a->d.allele[1],b->d.allele[1],b1)==0);
        }
    }

    bool bcf1_all_equal(bcf1_t *a, bcf1_t *b)
    {
        bcf_unpack(a, BCF_UN_ALL);
        bcf_unpack(b, BCF_UN_ALL);
        if (a == nullptr || b == nullptr)
        {
            die(" (bcf1_equal: tried to compare NULL bcf1_t");
        }
        if (a->rid != b->rid)
        {
            return (false);
        }
        else if (a->pos != b->pos)
        {
            return (false);
        }
        else
        {
            if(a->n_allele!=b->n_allele)
                return(false);

            for (size_t i = 0; i < a->n_allele; i++)
            {
                if (strcmp(a->d.allele[i], b->d.allele[i]))
                {
                    return (false);
                }
            }
        }
        return (true);
    }


    bool bcf1_less_than(bcf1_t *a, bcf1_t *b)
    {

        if (a == NULL || b == NULL)
        {
            die(" (bcf1_less_than: tried to compare NULL bcf1_t");
        }

        if (a->rid < b->rid)
        {
            return (true);
        }
        if (a->rid > b->rid)
        {
            return (false);
        }

        if (a->pos < b->pos)
        {
            return (true);
        }

        if (a->pos == b->pos)
        {
            if (get_variant_rank(a) == get_variant_rank(b))
            {
                size_t a1,b1,a2,b2;
                right_trim(a->d.allele[0],a->d.allele[1],a1,b1);
                right_trim(b->d.allele[0],b->d.allele[1],a2,b2);
                if(a1>a2 || b1>b2) return(false);
                int dif1 = strncmp(a->d.allele[0],b->d.allele[0],a1);
                int dif2 = strncmp(a->d.allele[1],b->d.allele[1],b1);
                if(dif1==0 && dif2==0) return(false);
                return(min(dif1,dif2)<0);
            }
            else
            {
                return (get_variant_rank(a) < get_variant_rank(b));
            }
        }
        return(false);
    }

    bool  bcf1_greater_than(bcf1_t *a, bcf1_t *b)
    {
        return (!bcf1_equal(a, b) && !bcf1_less_than(a, b));
    }

    bool  bcf1_leq(bcf1_t *a, bcf1_t *b)
    {
        return (!(bcf1_greater_than(a, b)));
    }

    bool bcf1_geq(bcf1_t *a, bcf1_t *b)
    {
        return (!(bcf1_less_than(a, b)));
    }

    bool bcf1_not_equal(bcf1_t *a, bcf1_t *b)
    {
        return (!(bcf1_equal(a, b)));
    }


    size_t get_number_of_likelihoods(int ploidy, int num_allele)
    {
        assert(ploidy == 1 || ploidy == 2);
        return (ploidy == 1 ? num_allele : (num_allele) * (1 + num_allele) / 2);
    }

    int phred(float l)
    {
        return ((int) (-10 * log10(l)));
    }

    float unphred(int pl)
    {
        return ((float) pow(10., -pl / 10.));
    }

    int factorial(int x)
    {
        int ret = 1;
        for (int i = 2; i <= x; i++)
        {
            ret *= i;
        }
        return (ret);
    }

    int choose(int n, int k)
    {
        if (k >= 0 && k <= n)
        {
            return (factorial(n) / (factorial(n - k) * factorial(k)));
        }
        else
        {
            return (0);
        }
    }

    //gets the index of a genotype likelihood for ploidy == 2
    int get_gl_index(int g0, int g1)
    {
        assert(g0 >= 0 && g1 >= 0);
        int a = g0 <= g1 ? g0 : g1;
        int b = g0 > g1 ? g0 : g1;
        return (choose(a, 1) + choose(b + 1, 2));
    }


    void print_variant(bcf_hdr_t const *header, bcf1_t *record)
    {
        bcf_unpack(record, BCF_UN_ALL);
        std::cerr << bcf_hdr_id2name(header, record->rid) << ":" << record->pos + 1 << ":" << record->d.allele[0];
        for (int i = 1; i < record->n_allele; i++)
        {
            std::cerr << ":" << record->d.allele[i];
        }
        int32_t *gt= nullptr,*ad= nullptr,*pl= nullptr;
        int num_gt=0,num_ad=0,num_pl=0;
        int ret = bcf_get_genotypes(header,record,&gt,&num_gt);
        if(ret>0)
        {
            std::cerr << "\tGT=";
            for (int i = 0; i < ret; i++)
            {
                if (i > 0) std::cerr << "/";
                std::cerr << bcf_gt_allele(gt[i]);
            }
        }

        ret=bcf_get_format_int32(header,record,"AD",&ad,&num_ad);
        if(ret>0)
        {
            std::cerr << "\tAD=";
            for (int i = 0; i < ret; i++)
            {
                if (i > 0) std::cerr << ",";
                std::cerr << ad[i];
            }
        }

        ret=bcf_get_format_int32(header,record,"PL",&pl,&num_pl);
        if(ret>0)
        {
            std::cerr << "\tPL=";
            for (int i = 0; i < ret; i++)
            {
                if (i > 0) std::cerr << ",";
                std::cerr << pl[i];
            }
        }
        free(gt);
        free(pl);
        free(ad);
        std::cerr << std::endl;
    }

    void print_variant(bcf1_t *record)
    {
        std::cerr << record->rid << ":" << record->pos + 1 << ":" << record->d.allele[0];
        for (int i = 1; i < record->n_allele; i++)
        {
            std::cerr << ":" << record->d.allele[i];
        }
        std::cerr << std::endl;
    }

    int get_ploidy(bcf_hdr_t *header, bcf1_t *record)
    {
        int *gt = NULL, ngt = 0;
        int ploidy = bcf_get_genotypes(header, record, &gt, &ngt);
        free(gt);
        return (ploidy);
    }

    int bcf1_get_one_info_float(bcf_hdr_t *header, bcf1_t *record, const char *tag,float & output)
    {
        bcf_unpack(record,BCF_UN_INFO);
        float *ptr=&output;
        int nval=1;
        int ret = bcf_get_info_float(header,record,tag,&ptr,&nval);
        if(ret>1)  die("bcf1_get_one_info_float:"+(string)tag+" more than one value returned");
        return(ret);
    }

    int bcf1_get_one_format_float(bcf_hdr_t *header, bcf1_t *record, const char *tag,float &output)
    {
        bcf_unpack(record, BCF_UN_FMT);
        float *ptr=&output;
        int nval=1;
        int ret =   bcf_get_format_float(header,record,tag,&ptr,&nval);
        if(ret>1)  die("bcf1_get_one_format_float:"+(string)tag+" more than one value returned");
        return(ret);
    }

    int bcf1_get_one_info_int(bcf_hdr_t *header, bcf1_t *record, const char *tag,int32_t & output)
    {
        bcf_unpack(record, BCF_UN_INFO);
        int32_t *ptr=&output;
        int nval=1;
        int ret = bcf_get_info_int32(header,record,tag,&ptr,&nval);
        if(ret>1)  die("bcf1_get_one_info_int: "+(string)tag+"more than one value returned");
        return(ret);
    }

    int bcf1_get_one_format_int(bcf_hdr_t *header, bcf1_t *record, const char *tag,int32_t &output)
    {
        bcf_unpack(record, BCF_UN_FMT);
        int32_t *ptr=&output;
        int nval=1;
        int ret = bcf_get_format_int32(header,record,tag,&ptr,&nval);
        if(ret>1)  die("bcf1_get_one_format_int:"+(string)tag+" more than one value returned");
        return(ret);
    }

    int bcf1_allele_swap(bcf_hdr_t *header, bcf1_t *record, int a,int b)
    {
        assert(a>0 && b>0);
        assert(a<record->n_allele && b<record->n_allele);
        bcf_unpack(record,BCF_UN_ALL);
        int32_t *format_ad=nullptr,*format_pl=nullptr,*gt= nullptr;
        int num_ad=0,num_pl=0,num_gt=0;
        vector<int> allele_map;//maps old alleles to new alleles
        for(int i=0;i<record->n_allele;i++) allele_map.push_back(i);
        allele_map[a]=b;
        allele_map[b]=a;

        //update the ALT
        char **new_alleles = (char **)malloc(record->n_allele*sizeof(char *));
        for(int i=0;i<record->n_allele;i++)
        {
            new_alleles[allele_map[i]] = (char *)malloc(strlen(record->d.allele[i])+1);
            strcpy(new_alleles[allele_map[i]],record->d.allele[i]);
        }
        bcf_update_alleles(header, record, (const char **) new_alleles, record->n_allele);

        //AD
        assert(bcf_get_format_int32(header,record,"AD",&format_ad,&num_ad)==record->n_allele);
        std::swap(format_ad[a],format_ad[b]);
        bcf_update_format_int32(header,record,"AD",format_ad,record->n_allele);


        //GT
        int ploidy=bcf_get_genotypes(header,record,&gt,&num_gt);
        assert(ploidy==1 || ploidy==2);
        bool phased = ploidy==1 ? bcf_gt_is_phased(gt[0]) : bcf_gt_is_phased(gt[0])&&bcf_gt_is_phased(gt[1]);
        for(int i=0;i<ploidy;i++)
        {
            if(bcf_gt_allele(gt[i])==a)
                gt[i] = phased ? bcf_gt_phased(b) : bcf_gt_unphased(b);
            else if(bcf_gt_allele(gt[i])==b)
                gt[i] = phased ? bcf_gt_phased(a) : bcf_gt_unphased(a);
        }
        bcf_update_genotypes(header,record,gt,ploidy);

        //PL
        assert(bcf_get_format_int32(header,record,"PL",&format_pl,&num_pl)==(int)ggutils::get_number_of_likelihoods(ploidy,record->n_allele));
        vector<int> tmp_pl(format_pl,format_pl+num_pl);
        for(int i=0;i<record->n_allele;i++)
            for(int j=i;j<record->n_allele;j++)
                format_pl[ggutils::get_gl_index(allele_map[i],allele_map[j])] = tmp_pl[ggutils::get_gl_index(i,j)];

        bcf_update_format_int32(header,record,"PL",format_pl,num_pl);

        //memory cleanup
        for(int i=0;i<record->n_allele;i++) free(new_alleles[i]);
        free(new_alleles);
        free(format_ad);
        free(format_pl);
        return(1);
    }
}
