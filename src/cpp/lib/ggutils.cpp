#include <htslib/vcf.h>
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

    int32_t *zeros(int n)
    {
        int32_t *ret = (int *) malloc(sizeof(int32_t) * n);
        std::fill(ret,ret+n,0);
        return (ret);
    }

    int32_t *assign_bcf_int32_missing(int32_t *ptr,size_t n)
    {
        ptr = (int32_t *)realloc(ptr,sizeof(int32_t) * n);
        std::fill(ptr,ptr+n,bcf_int32_missing);
        return(ptr);
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
        return bcf_is_snp(record);
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


    bool bcf1_less_than(bcf1_t *a, bcf1_t *b)
    {
        assert(a->n_allele>1 && b->n_allele>1);
        assert(a!=nullptr && b!=nullptr);

        if (a->rid < b->rid)
            return(true);
        if (a->rid > b->rid)
            return(false);
        if (a->pos < b->pos)
            return(true);
        if (a->pos > b->pos)
            return(false);
        if(get_variant_rank(a) < get_variant_rank(b))
            return(true);
        if(get_variant_rank(a) > get_variant_rank(b))
            return(false);
        if(bcf1_equal(a,b))
            return(false);
        size_t a1,b1,a2,b2;
        right_trim(a->d.allele[0],a->d.allele[1],a1,b1);
        right_trim(b->d.allele[0],b->d.allele[1],a2,b2);
        if(a1<a2)
            return(true);
        if(b1<b2)
            return(true);
        if(strncmp(a->d.allele[0],b->d.allele[0],min(a1,a2))<0)
            return(true);
        if(strncmp(a->d.allele[1],b->d.allele[1],min(b1,b2))<0)
            return(true);
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

    //FIXME: One needs to be careful of overlow here. I have just used a recursive implementation. There is probably a better way to do this (check Rmath.h for example).
    int choose(int n, int k)
    {
        if(n==0)
            return(0);
        if(k==0 || k==n)
            return(1);
        else
            return(choose(n-1,k-1)+choose(n-1,k));
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
	std::cerr<<record2string(header,record) <<std::endl;
    }
	
    string record2string(bcf_hdr_t const *header, bcf1_t *record)
    {
        bcf_unpack(record, BCF_UN_ALL);
	std::stringstream ss;
        ss << bcf_hdr_id2name(header, record->rid) << ":" << record->pos + 1 << ":" << record->d.allele[0];
        for (int i = 1; i < record->n_allele; i++)
        {
            ss << ":" << record->d.allele[i];
        }
	return ss.str();
        // int32_t *gt= nullptr,*ad= nullptr,*pl= nullptr;
        // int num_gt=0,num_ad=0,num_pl=0;
        // int ret = bcf_get_genotypes(header,record,&gt,&num_gt);
        // if(ret>0)
        // {
        //     ss << "\tGT=";
        //     for (int i = 0; i < ret; i++)
        //     {
        //         if (i > 0) ss << "/";
        //         ss << bcf_gt_allele(gt[i]);
        //     }
        // }

        // ret=bcf_get_format_int32(header,record,"AD",&ad,&num_ad);
        // if(ret>0)
        // {
        //     ss << "\tAD=";
        //     for (int i = 0; i < ret; i++)
        //     {
        //         if (i > 0) ss << ",";
        //         ss << ad[i];
        //     }
        // }

        // ret=bcf_get_format_int32(header,record,"PL",&pl,&num_pl);
        // if(ret>0)
        // {
        //     ss << "\tPL=";
        //     for (int i = 0; i < ret; i++)
        //     {
        //         if (i > 0) ss << ",";
        //         ss << pl[i];
        //     }
        // }
        // free(gt);
        // free(pl);
        // free(ad);
        // ss << std::endl;
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

    int bcf1_get_one_info_float(const bcf_hdr_t *header, bcf1_t *record, const char *tag,float & output)
    {
        bcf_unpack(record,BCF_UN_INFO);
        float *ptr=&output;
        int nval=1;
        int ret = bcf_get_info_float(header,record,tag,&ptr,&nval);
        if(ret>1) die("bcf1_get_one_info_float:"+(string)tag+" more than one value returned");
        if(ret<0) output=bcf_float_missing;
        return(ret);
    }

    int bcf1_get_one_format_float(const bcf_hdr_t *header, bcf1_t *record, const char *tag,float &output)
    {
        bcf_unpack(record, BCF_UN_FMT);
        float *ptr=&output;
        int nval=1;
        int ret =   bcf_get_format_float(header,record,tag,&ptr,&nval);
        if(ret>1) die("bcf1_get_one_format_float:"+(string)tag+" more than one value returned");
        if(ret<0) output=bcf_float_missing;
        return(ret);
    }

    int bcf1_get_one_info_int(const bcf_hdr_t *header, bcf1_t *record, const char *tag,int32_t & output)
    {
        bcf_unpack(record, BCF_UN_INFO);
        int32_t *ptr=&output;
        int nval=1;
        int ret = bcf_get_info_int32(header,record,tag,&ptr,&nval);
        if(ret>1) die("bcf1_get_one_info_int: "+(string)tag+"more than one value returned");
        if(ret<0) output=bcf_int32_missing;
        return(ret);
    }

    int bcf1_get_one_format_int(const bcf_hdr_t *header, bcf1_t *record, const char *tag,int32_t &output)
    {
        bcf_unpack(record, BCF_UN_FMT);
        int32_t *ptr=&output;
        int nval=1;
        int ret = bcf_get_format_int32(header,record,tag,&ptr,&nval);
        if(ret>1) die("bcf1_get_one_format_int:"+(string)tag+" more than one value returned");
        if(ret<0) output=bcf_int32_missing;
        return(ret);
    }

    int bcf1_allele_swap(bcf_hdr_t *header, bcf1_t *record, int a,int b)
    {
        assert(a>0 && b>0);
        assert(a<record->n_allele && b<record->n_allele);
        bcf_unpack(record,BCF_UN_ALL);
        int32_t *format_ad=nullptr,*format_adr=nullptr,*format_adf=nullptr,*format_pl=nullptr,*gt= nullptr;
        int num_adr=0,num_adf=0,num_ad=0,num_pl=0,num_gt=0;
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

        if(bcf_get_format_int32(header,record,"ADF",&format_adf,&num_adf)==record->n_allele)
        {
            std::swap(format_adf[a], format_adf[b]);
            bcf_update_format_int32(header, record, "ADF", format_adf, record->n_allele);
        }

        if(bcf_get_format_int32(header,record,"ADR",&format_adr,&num_adr)==record->n_allele)
        {
            std::swap(format_adr[a], format_adr[b]);
            bcf_update_format_int32(header, record, "ADR", format_adr, record->n_allele);
        }

        //GT
        int ploidy=bcf_get_genotypes(header,record,&gt,&num_gt);
	if(ploidy==2 && gt[1]==bcf_int32_vector_end) ploidy=1;
	    
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
        int status = bcf_get_format_int32(header,record,"PL",&format_pl,&num_pl);
        if(status>0)
        {
            if(status!=(int)ggutils::get_number_of_likelihoods(ploidy,record->n_allele))
	    {
		ggutils::print_variant(header,record);
		ggutils::die("problem with sample "+(string)header->samples[0]);
	    }
            vector<int> tmp_pl(format_pl,format_pl+num_pl);
            for(int i=0;i<record->n_allele;i++)
            {
                if(ploidy==1)
                    format_pl[allele_map[i]] = tmp_pl[i];
                else
                    for(int j=i;j<record->n_allele;j++)
                    {
                        format_pl[ggutils::get_gl_index(allele_map[i], allele_map[j])] = tmp_pl[ggutils::get_gl_index(i, j)];
                    }
            }
            bcf_update_format_int32(header,record,"PL",format_pl,num_pl);
        }

        //memory cleanup
        for(int i=0;i<record->n_allele;i++) free(new_alleles[i]);
        free(new_alleles);
        free(format_ad);
        free(format_adr);
        free(format_adf);
        free(format_pl);
        free(gt);
        return(1);
    }

    vcf_data_t::vcf_data_t(size_t ploidy, size_t num_allele, size_t num_sample)
    {
        this->ploidy=ploidy;
        this->num_allele=num_allele;
        this->num_sample=num_sample;

        gq=(int32_t*)malloc(num_sample*sizeof(int32_t));
        gqx=(int32_t*)malloc(num_sample*sizeof(int32_t));
        dp=(int32_t*)malloc(num_sample*sizeof(int32_t));
        dpf=(int32_t*)malloc(num_sample*sizeof(int32_t));
        ps=(int32_t*)malloc(num_sample*sizeof(int32_t));

        num_ad=num_allele*num_sample;
        ad = (int32_t *)malloc(num_ad*sizeof(int32_t));
        adf = (int32_t *)malloc(num_ad*sizeof(int32_t));
        adr = (int32_t *)malloc(num_ad*sizeof(int32_t));

        num_pl = get_number_of_likelihoods(ploidy,num_allele)*num_sample;
        pl = (int32_t *)malloc(num_pl*sizeof(int32_t));
        gt = (int32_t *)malloc(num_sample*ploidy*sizeof(int32_t));
    }

    void vcf_data_t::resize(size_t num_alleles)
    {
        if(num_alleles!=this->num_allele)
        {
            this->num_allele=num_alleles;
            num_pl = ggutils::get_number_of_likelihoods(ploidy,num_allele)* num_sample;
            pl = (int32_t *) realloc(pl, num_pl * sizeof(int32_t));
            num_ad=num_allele*num_sample;
            ad = (int32_t *)realloc(ad,num_ad*sizeof(int32_t));
            adr = (int32_t *)realloc(adr,num_ad*sizeof(int32_t));
            adf = (int32_t *)realloc(adf,num_ad*sizeof(int32_t));
        }
    }

    void vcf_data_t::set_missing()
    {
        std::fill(gq,gq+num_sample,bcf_int32_missing);
        std::fill(gqx,gqx+num_sample,bcf_int32_missing);
        std::fill(dp,dp+num_sample,bcf_int32_missing);
        std::fill(dpf,dpf+num_sample,bcf_int32_missing);
        std::fill(ps,ps+num_sample,bcf_int32_missing);
        std::fill(gt,gt+ploidy*num_sample,bcf_gt_missing);
        std::fill(pl,pl+num_pl,bcf_int32_missing);
        std::fill(ad,ad+num_ad,bcf_int32_missing);
        std::fill(adf,adf+num_ad,bcf_int32_missing);
        std::fill(adr,adr+num_ad,bcf_int32_missing);
    }

    vcf_data_t::~vcf_data_t()
    {
        free(ad);
        free(adf);
        free(adr);
        free(pl);
        free(gt);
        free(gq);
        free(gqx);
        free(dp);
        free(dpf);
        free(ps);
    }

    int find_allele(bcf1_t *target,bcf1_t *query,int index)
    {
        bcf_unpack(query,BCF_UN_ALL);
        bcf_unpack(target,BCF_UN_ALL);
        assert(target!=nullptr);
        assert(query!=nullptr);
        assert(index>0 && index<query->n_allele);
        size_t rlen,alen;
        char *q_ref=query->d.allele[0];
        char *q_alt=query->d.allele[index];
        ggutils::right_trim(q_ref,q_alt,rlen,alen);
        char *t_ref=target->d.allele[0];
        for(int i=1;i<target->n_allele;i++)
        {
            size_t target_rlen,target_alen;
            char *t_alt=target->d.allele[i];
            ggutils::right_trim(t_ref,t_alt,target_rlen,target_alen);
            if(rlen==target_rlen && alen==target_alen)
                if(strncmp(q_ref,t_ref,rlen)==0 && strncmp(q_alt,t_alt,alen)==0)
                    return(i);
        }
        return(-1);
    }

    int add_allele(bcf_hdr_t *hdr,bcf1_t *dst,bcf1_t *src,int index)
    {
        bcf_unpack(dst,BCF_UN_ALL);
        bcf_unpack(src,BCF_UN_ALL);
        assert(index>0 && index<src->n_allele);
        int dst_index = ggutils::find_allele(dst,src,index);
        if(dst_index>0) 
            return(dst_index);
        size_t num_alleles = dst->n_allele+1;
        char **new_alleles = (char **)malloc(sizeof(char*) * num_alleles);
        size_t old_ref_len = strlen(dst->d.allele[0]);
        size_t q_rlen,q_alen;
        char *q_ref=src->d.allele[0];
        char *q_alt=src->d.allele[index];
        ggutils::right_trim(q_ref,q_alt,q_rlen,q_alen);
        if(old_ref_len==q_rlen)
        {
            for(size_t i=0;i<num_alleles-1;i++)
                new_alleles[i]=strdup(dst->d.allele[i]);
            new_alleles[num_alleles-1]=(char *)malloc(q_alen+1);
            strncpy(new_alleles[num_alleles-1],src->d.allele[index],q_alen);
            new_alleles[num_alleles-1][q_alen]='\0';
        }
        else if(old_ref_len>q_rlen)
        {
            for(size_t i=0;i<num_alleles-1;i++)
                new_alleles[i]=strdup(dst->d.allele[i]);
            size_t rightpad = old_ref_len - q_rlen;
            new_alleles[num_alleles-1] = (char *)malloc(q_alen + rightpad + 1);
            memcpy(new_alleles[num_alleles-1],q_alt,q_alen);
            memcpy(new_alleles[num_alleles-1]+q_alen,new_alleles[0]+old_ref_len-rightpad,rightpad);
            new_alleles[num_alleles-1][q_alen+rightpad]='\0';
        }
        else
        {
            size_t rightpad = q_rlen - old_ref_len;
            for(size_t i=0;i<num_alleles-1;i++)
            {
                int alen=strlen(dst->d.allele[i]);
                new_alleles[i]=(char *)malloc(alen+rightpad+1);
                memcpy(new_alleles[i],dst->d.allele[i],alen);
                memcpy(new_alleles[i]+alen,q_ref+old_ref_len,rightpad);
                new_alleles[i][alen+rightpad]='\0';
            }
            new_alleles[num_alleles-1] = (char *)malloc(q_alen+1);
            memcpy(new_alleles[num_alleles-1],q_alt,q_alen);
            new_alleles[num_alleles-1][q_alen]='\0';
        }
        bcf_update_alleles(hdr,dst,(const char**)new_alleles,num_alleles);
        for(size_t i=0;i<num_alleles;i++) free(new_alleles[i]);
        free(new_alleles);
        return(num_alleles-1);
    }

    void collapse_haploid_gls(std::vector< std::vector<int> > & pls,std::vector<int> & output)
    {
        output[0]=0;
        for(auto it=pls.begin();it!=pls.end();it++) output[0] += (*it)[0];
        for(auto it=pls.begin();it!=pls.end();it++)
        {
            for(size_t i=1;i<it->size();i++)
            {
                if((*it)[i]!=bcf_int32_missing)
                    output[i]=(*it)[i];
                else
                    output[i]+=(*it)[0];
            }
        }
    }

    void collapse_gls(int ploidy,int num_alleles,std::vector< std::vector<int> > & pls,std::vector<int> & output)
    {
        assert(ploidy==1 || ploidy==2);
        int num_sets = pls.size();
        assert(num_sets>0);
        size_t num_gls = pls[0].size();
        for(auto it=pls.begin();it!=pls.end();it++) assert(it->size()==num_gls);
        output.assign(ggutils::get_number_of_likelihoods(ploidy,num_alleles),bcf_int32_missing);
        if(ploidy==1)
        {
            collapse_haploid_gls(pls,output);
            return;
        }

        int g00=0;
        for(auto it=pls.begin();it!=pls.end();it++)
            g00 += (*it)[0];

        for(auto it1=pls.begin();it1!=pls.end();it1++)
        {
            int gi0 = 0;
            for(auto it2=pls.begin();it2!=pls.end();it2++)
                if(it1!=it2)
                    gi0 += (*it2)[0];
            
            for(int i=0;i<num_alleles;i++)
            {
                for(int j=i;j<num_alleles;j++)
                {
                    int index = ggutils::get_gl_index(i,j);
                    if(i==0 && j==0)
                    {
                        output[0]=g00;  // collapse all the GT=0/0 GLs
                    }
                    else if(i==j || i==0)
                    {
                        int gl = (*it1)[index];
                        if(gl != bcf_int32_missing)
                            output[index] = gl + gi0; //GT=i/i or GT=0/i
                    }
                    else //GT=i/j where i!=j && i>0 && j>0
                    {
                        int gij = 0;
                        for(auto it2=pls.begin();it2!=pls.end();it2++)
                        {
                            if((*it2)[index] == bcf_int32_missing)
                            {                             
                                int gi = (*it2)[ggutils::get_gl_index(0,i)];
                                int gj = (*it2)[ggutils::get_gl_index(0,j)];
                                if(gi != bcf_int32_missing) {gij += gi;}
                                else if(gj != bcf_int32_missing) {gij += gj;}
                                else {gij += (*it2)[0];}
                            }
                            else
                            {
                                gij += (*it2)[index];
                            }
                        }
                        output[index] = gij;
                    }
                }
            }
        }
         int min_pl = *std::min_element(output.begin(), output.end());
         for(size_t i=0;i<output.size();i++) output[i] -= min_pl;
    }

    std::string string_time()
    {
        std::time_t result = std::time(nullptr);
        char buffer[16];
        std::strftime(buffer,120,"%Y%m%d_%H%M%S",std::localtime(&result));
        return((std::string)buffer);
    }


    float median(int *x, int n)
    {
        assert(n>0);
        std::vector<int> work(x,x+n);
        return(inplace_median(work));
    }
    
    //TODO: I think a better way to do this is the quickselect algorithm: https://en.wikipedia.org/wiki/Quickselect
    float inplace_median(std::vector<int> & work)
    {
        size_t n = work.size();
        assert(n>0);
        if(n%2==1)
        {
            std::nth_element(work.begin(),work.begin()+n/2,work.end());
            return work[n/2];
        }
        else
        {
            std::nth_element(work.begin(),work.begin()+n/2,work.end());
            float a=work[n/2];
            std::nth_element(work.begin(),work.begin()+n/2-1,work.end());
            a+=work[n/2-1];
            return a/2.;
        }        
    }   
    
    void fisher_sb_test(int *adf,int *adr,int num_allele,std::vector<float> & output,float maxret)
    {
        assert(num_allele>1);
        output.resize(num_allele-1);
        double left,right,two_sided;
        for(int i=1;i<num_allele;i++)
        {
            kt_fisher_exact(adr[0],adr[i],adf[0],adf[i],&left,&right,&two_sided);
            output[i-1] = -log10(two_sided);
	        if(output[i-1]== -0.0) output[i-1]=0.;//gets rid of silly -0.0
            if(output[i-1]>maxret) output[i-1]=maxret;
        }
    }
    
    bool is_valid_strelka_record(bcf_hdr_t const *header, bcf1_t *record)
    {
	int32_t *int_ptr=nullptr;
	int num_int=0;
	int status = bcf_get_format_int32(header, record, "AD", &int_ptr, &num_int);
	if(status>0) free(int_ptr);
	return(status == record->n_allele);
    }
}
