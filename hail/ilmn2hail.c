/*  ilmn2hail.c: converts haploid calls to diploid calls by padding with reference
                 genotypes. If FORMAT/PL is present, dummy 0,255,255 or 255,0,255 
        		 values will be inserted accordingly (for example).
                 Sets arrays INFO/ADF and INFO/ADR to missing iff one of their entries
                 is missing.

    Copyright (C) 2017 Illumina

    Authors: Ole Schulz-Trieglaff <oschulz-trie@illumina.com>, Jared O'Connell <jared.oconnell@gmail.com>

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
    THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.  */

#include <stdio.h>
#include <stdlib.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <math.h>

bcf_hdr_t *in_hdr, *out_hdr;
int *gt = NULL,*gt_out=NULL, ngt = 0,nsample;
int *adf = NULL,*adf_out=NULL, nadf = 0;
int *adr = NULL,*adr_out=NULL, nadr = 0;
int32_t *pl=NULL,*pl_out,npl=0;

const char *about(void)
{
    return "converts haploid calls into diploid calls by padding with reference genotypes\n";
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    in_hdr  = in;
    out_hdr = out;
    nsample =  bcf_hdr_nsamples(in_hdr);
    gt_out = (int *)malloc(sizeof(int32_t)*2*nsample);
    adf_out = (int *)malloc(sizeof(int32_t)*2*nsample);
    adr_out = (int *)malloc(sizeof(int32_t)*2*nsample);
    return 0;
}

int get_number_of_likelihoods(int ploidy, int num_allele)
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

bcf1_t *process(bcf1_t *rec)
{
    int i,j;
    int n_pl  =  bcf_get_format_int32(in_hdr,rec,"PL",&pl,&npl);
    if (n_pl < 0) {
        fprintf(stderr,"Read %d genotype likelihoods (FORMAT/PL) but expected at least one. Aborting\n",n_pl);
        assert(0);
    }
    int n_adf =  bcf_get_format_int32(in_hdr,rec,"ADF",&adf,&nadf);
    if (n_adf < 0) {
        fprintf(stderr,"Error reading allelic depth (FORMAT/ADF). Aborting\n");
        assert(0);
    }
    int n_adr =  bcf_get_format_int32(in_hdr,rec,"ADR",&adr,&nadr);
    if (n_adr < 0) {
        fprintf(stderr,"Error reading allelic depth (FORMAT/ADR). Aborting\n");
        assert(0);
    }


    int gt_ret = bcf_get_genotypes(in_hdr,rec,&gt,&ngt);
    if (gt_ret!=2*nsample && gt_ret!=nsample) {
        fprintf(stderr,"Read %d genotypes but expected %d or %d. Aborting\n",gt_ret,(2*nsample),nsample);
        assert(0);
    }
    int ploidy=2;
    int nal=rec->n_allele;
    int num_pl_per_sample = get_number_of_likelihoods(ploidy,nal);
    pl_out = (int32_t *)realloc(pl_out,nsample*num_pl_per_sample*sizeof(int32_t*));
				
    for(i=0;i<nsample;i++) 
    {
        for(j=0;j<num_pl_per_sample;j++)
            pl_out[i*num_pl_per_sample + j] = 255;

        if(gt_ret==nsample)//all haploid
        {
            if(bcf_gt_is_missing(gt[i]))
            {
            gt_out[2*i]=bcf_gt_missing;
            gt_out[2*i+1]=bcf_gt_missing;
            }
            else
            {
            gt_out[2*i] = bcf_gt_unphased(0);
            gt_out[2*i+1] = bcf_gt_unphased(bcf_gt_allele(gt[i]));
            pl_out[i*num_pl_per_sample + get_gl_index(bcf_gt_allele(gt_out[2*i]),bcf_gt_allele(gt_out[2*i+1]))] = 0 ;		
            }
        }
        else if(gt[2*i+1]==bcf_int32_vector_end)//this sample is haploid
        {
            if(bcf_gt_is_missing(gt[2*i]))
            {
            gt_out[2*i]=bcf_gt_missing;
            gt_out[2*i+1]=bcf_gt_missing;
            }
            else
            {		
            gt_out[2*i] = bcf_gt_unphased(0);
            gt_out[2*i+1] = bcf_gt_unphased(bcf_gt_allele(gt[2*i]));
            pl_out[i*num_pl_per_sample + get_gl_index(bcf_gt_allele(gt_out[2*i]),bcf_gt_allele(gt_out[2*i+1]))] = 0 ;
            }
        }
        else
        {
            gt_out[2*i] = gt[2*i];
            gt_out[2*i+1] = gt[2*i+1];
            if(!bcf_gt_is_missing(gt[2*i]))
            {
            for(j=0;j<num_pl_per_sample;j++)
                pl_out[i*num_pl_per_sample + j] = pl[i*num_pl_per_sample + j];
            }
        }
    
        // If one of the entries of FORMAT/ADR or FORMAT/ADR is missing,
        // set the whole array to missing
        //http://discuss.hail.is/t/matrix-table-error/485/2
        if (adf[2*i]==bcf_int32_missing || adf[2*i+1]==bcf_int32_missing) {
            adf_out[2*i]   = bcf_int32_missing;
            adf_out[2*i+1] = bcf_int32_vector_end;
        } else {
            adf_out[2*i]   = adf[2*i];
            adf_out[2*i+1] = adf[2*i+1];

        }
        if (adr[2*i]==bcf_int32_missing || adr[2*i+1]==bcf_int32_missing) {
            adr_out[2*i]   = bcf_int32_missing;
            adr_out[2*i+1] = bcf_int32_vector_end;
        } else {
            adr_out[2*i]   = adr[2*i];
            adr_out[2*i+1] = adr[2*i+1];
        }
        
    }
    
    bcf_update_genotypes(out_hdr,rec,gt_out,2*nsample);
    bcf_update_format_int32(out_hdr,rec,"PL",pl_out,nsample*num_pl_per_sample);
    bcf_update_format_int32(out_hdr,rec,"ADF",adf_out,2*nsample);
    bcf_update_format_int32(out_hdr,rec,"ADR",adr_out,2*nsample);

    return rec;
}

void destroy(void)
{
    free(gt);
    free(pl);
    free(gt_out);
    free(pl_out);
}
