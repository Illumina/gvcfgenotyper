#include "utils.hh"

using namespace std;

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
    for (int i = 1; i < input.size(); i++)
    {
        ret += ",";
        ret += input[i];
    }
    return (ret);
}

vector<int> match(const vector<string> &x, const vector<string> &y)
{
    vector<int> ret;
    for (int i = 0; i < y.size(); i++)
    {
        int j = 0;
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
    kstring_t splitme = {0,0,nullptr};
    bcf_hdr_format(src,1,&splitme);
    strsplit(splitme.s, '\n', parseme);
    string contigs = "";
    for (int i = 0; i < parseme.size(); i++)
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
        die("problem opening" + ((string) fname));
    }
    string tmp;
    output.clear();
    while (getline(ifile,tmp))
    {
        output.push_back(tmp);
    }

    return (output.size());
}


bool is_snp(bcf1_t *record)
{
    assert(record->n_allele>1);
    bcf_unpack(record,BCF_UN_ALL);
    return(bcf_get_variant_type(record,1)&VCF_SNP);
}

bool is_deletion(bcf1_t *record)
{
    bcf_unpack(record,BCF_UN_ALL);
    int l1 = strlen(record->d.allele[0]);
    int l2 = strlen(record->d.allele[1]);
    return (bcf_get_variant_type(record,1) == VCF_INDEL && l2<l1);
}

bool is_insertion(bcf1_t *record)
{
    bcf_unpack(record,BCF_UN_ALL);
    int l1 = strlen(record->d.allele[0]);
    int l2 = strlen(record->d.allele[1]);
    return (bcf_get_variant_type(record,1) == VCF_INDEL && l2>l1);
}

bool is_complex(bcf1_t *record)
{
    return(!is_snp(record) && !is_deletion(record) && !is_insertion(record));
}

int get_variant_rank(bcf1_t *record)
{
    if(is_complex(record)||is_snp(record))
    {
        return(0);
    }
    if(is_insertion(record))
    {
        return(1);
    }
    if(is_deletion(record))
    {
        return(2);
    }
    die("bad variant");
    return(-1);
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

bool bcf1_equal(bcf1_t *a, bcf1_t *b)
{
    bcf_unpack(a,BCF_UN_ALL);
    bcf_unpack(b,BCF_UN_ALL);
    if (a == nullptr || b == nullptr)
    {
        die("bcf1_equal: tried to compare NULL bcf1_t");
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
        for (int i = 0; i < min(a->n_allele,b->n_allele); i++)
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
        die("bcf1_less_than: tried to compare NULL bcf1_t");
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
        if(get_variant_rank(a)==get_variant_rank(b))
        {
            for (int i = 0; i < min(a->n_allele, b->n_allele); i++)
            {
                int val = strcmp(a->d.allele[i], b->d.allele[i]);
                if (val < 0)
                {
                    return (true);
                }
                if (val > 0)
                {
                    return (false);
                }
            }
        }
        else
        {
            return(get_variant_rank(a)<get_variant_rank(b));
        }
    }
    return (false);
}

bool bcf1_greater_than(bcf1_t *a, bcf1_t *b)
{
    return (!bcf1_equal(a, b) && !bcf1_less_than(a, b));
}

bool bcf1_leq(bcf1_t *a, bcf1_t *b)
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
    return (ploidy == 1 ? ploidy : (num_allele) * (1 + num_allele) / 2);
}

int  phred(float l)
{
    return ((int) (-10 * log10(l)));
}

float  unphred(int pl)
{
    return ((float) pow(10., -pl / 10.));
}

int  factorial(int x)
{
    int ret = 1;
    for (int i = 2; i <= x; i++)
    {
        ret *= i;
    }
    return (ret);
}

int  choose(int n, int k)
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
int  get_gl_index(int g0, int g1)
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
