#pragma once

#define __STDC_LIMIT_MACROS

#include <stdint.h>
#include <limits>

#include "math.h"
#include <set>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdlib.h>

using namespace std;

extern "C" {
#include <htslib/bgzf.h>
#include <htslib/tbx.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <htslib/synced_bcf_reader.h>
}


static int *zeros(int n)
{
    int *ret = (int *)malloc(sizeof(int)*n);
    for(int i=0;i<n;i++)
    {
        ret[i]=0;
    }
    return(ret);
}

template<typename T>
void assign(int n,T val,T *x)
{
    for(size_t i=0;i<n;i++)
    {
        x[i]=val;
    }
};

bool fileexists(const string &fname);

template<typename T>
void newMatrix(int nrow, int ncol, T **&out)
{
    out = new T *[nrow];
    for (int i = 0; i < nrow; i++)
        out[i] = new T[ncol];
    //  cout << out[nrow-1][ncol-1]<<endl;
};


template<typename T>
void printMatrix(int nrow, int ncol, T *&out)
{
    for (int i = 0; i < nrow; i++)
    {
        for (int j = 0; j < ncol; j++)
        {
            cout << out[i * ncol + j] << " ";
        }
        cout << endl;
    }
    cout << endl;
};

template<typename T>
void delMatrix(int nrow, int ncol, T **&out)
{
    for (int i = 0; i < nrow; i++)
    {
        delete[] out[i];
    }
    delete[] out;
};


//retruns first location where find_this occurs in_that (-1 if not found);
template<typename T>
int find(T find_this, vector<T> in_that)
{
    for (int i = 0; i < in_that.size(); i++)
    {
        if (in_that[i] == find_this)
        {
            return (i);
        }
    }
    return (-1);
};


inline float argmax(float *P, int K, uint &maxi, uint &maxj)
{
    //  assert(!isnan(*ptr));
    float *ptr = P;
    float maxval = *ptr;
    maxi = 0;
    maxj = 0;
    for (int i = 0; i < K; i++)
    {
        for (int j = 0; j < K; j++)
        {
            if (*ptr > maxval)
            {
                maxval = *ptr;
                maxi = i;
                maxj = j;
            }
            ptr++;
        }
    }
    return (maxval);
}

string percent(int num, int den);

inline void die(const string &s)
{
    cerr << "ERROR: " << s << "\nExiting..." << endl;
    exit(1);
}

int strsplit(const string &input, const char split, vector<string> &out);

string join(const vector<string> &input, const string &delim);

vector<int> match(const vector<string> &x, const vector<string> &y);

inline void warn(const string &s)
{
    cerr << "WARNING: " << s << endl;
}

int copy_contigs(const bcf_hdr_t *src, bcf_hdr_t *dst);

//simple dumps text from fname into output
int read_text_file(const string &fname, vector<string> &output);


static bool is_snp(bcf1_t *record)
{
    return (record->n_allele == 2 && strlen(record->d.allele[0]) == 1 && strlen(record->d.allele[1]) == 1);
}

static int get_end_of_gvcf_block(bcf_hdr_t *header, bcf1_t *record)
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

static int get_end_of_variant(bcf1_t *record)
{
    return (record->pos + strlen(record->d.allele[0]) - 1);
}

static bool bcf1_equal(const bcf1_t *a, const bcf1_t *b)
{
    if (a == NULL || b == NULL)
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
    else if (a->n_allele != b->n_allele)
    {
        return (false);
    }
    else
    {
        for (int i = 0; i < a->n_allele; i++)
        {
            if (strcmp(a->d.allele[i], b->d.allele[i]))
            {
                return (false);
            }
        }
    }
    return (true);
}


static bool bcf1_less_than(const bcf1_t *a, const bcf1_t *b)
{
    if (a == NULL || b == NULL)
    {
        die("bcf1_less_than: tried to compare NULL bcf1_t");
    }

    if (a->rid < b->rid)
    {
        return (true);
    }

    if (a->pos < b->pos)
    {
        return (true);
    }

    if (a->pos == b->pos)
    {
        if (a->n_allele < b->n_allele)
        {
            return (true);
        }
        else
        {
            for (int i = 0; i < a->n_allele; i++)
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
    }
    return (false);
}

static bool bcf1_greater_than(const bcf1_t *a, const bcf1_t *b)
{
    return (!bcf1_equal(a, b) && !bcf1_less_than(a, b));
}

static bool bcf1_leq(const bcf1_t *a, const bcf1_t *b)
{
    return (!(bcf1_greater_than(a, b)));
}

static bool bcf1_geq(const bcf1_t *a, const bcf1_t *b)
{
    return (!(bcf1_less_than(a, b)));
}

static bool bcf1_not_equal(const bcf1_t *a, const bcf1_t *b)
{
    return (!(bcf1_equal(a, b)));
}

static void print_variant(bcf_hdr_t *header, bcf1_t *record)
{
    bcf_unpack(record, BCF_UN_ALL);
    std::cerr << bcf_hdr_id2name(header, record->rid) << ":" << record->pos + 1 << ":" << record->d.allele[0];
    for (int i = 1; i < record->n_allele; i++)
    {
        std::cerr << ":" << record->d.allele[i];
    }
    std::cerr << std::endl;
}

static size_t get_number_of_likelihoods(int ploidy,int num_allele)
{
    assert(ploidy==1 || ploidy==2);
    return(ploidy==1 ? ploidy : (num_allele) * (1 + num_allele) / 2);
}

static int inline phred(float l)
{
    return((int)(-10*log10(l)));
}

static float inline unphred(int pl)
{
    return((float)pow(10.,-pl/10.));
}

static int inline factorial(int x)
{
    int ret = 1;
    for(int i=2;i<=x;i++)
    {
        ret *= i;
    }
    return(ret);
}

static int inline choose(int n,int k)
{
    if(k>=0 && k<=n)
    {
        return (factorial(n) / (factorial(n - k) * factorial(k)));
    }
    else
    {
        return(0);
    }
}

//gets the index of a genotype likelihood for ploidy == 2
static int inline get_gl_index(int g0,int g1)
{
    assert(g0>=0 && g1>=0);
    int a = g0<=g1 ? g0 : g1;
    int b = g0>g1 ? g0 : g1;
    return(choose(a,1) + choose(b + 1,2));
}
