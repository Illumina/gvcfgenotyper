//This is a modified version of vcfnorm.c from bcftools. The original copyright notice is below. Jared O'Connell <joconnell@illumina>

/*  vcfnorm.c -- Left-align and normalize indels.

    Copyright (C) 2013-2014 Genome Research Ltd.

    Author: Petr Danecek <pd3@sanger.ac.uk>

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
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    THE SOFTWARE.  */

#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/faidx.h>
#include "vcfnorm.h"

#define CHECK_REF_EXIT 0
#define CHECK_REF_WARN 1
#define CHECK_REF_SKIP 2
#define CHECK_REF_FIX  4

#define MROWS_SPLIT 1
#define MROWS_MERGE  2


static inline int replace_iupac_codes(char *seq, int nseq)
{
    // Replace ambiguity codes with N for now, it awaits to be seen what the VCF spec codifies in the end
    int i, n = 0;
    for (i = 0; i < nseq; i++)
    {
        seq[i] = toupper(seq[i]);
        char c = seq[i];
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T')
        {
            seq[i] = 'N';
            n++;
        }
    }
    return n;
}

#define ERR_DUP_ALLELE      -2
#define ERR_REF_MISMATCH    -1
#define ERR_OK              0
#define ERR_SYMBOLIC        1

int realign(args_t *args, bcf1_t *line,bcf_hdr_t *hdr)
{
    bcf_unpack(line, BCF_UN_STR);

    // Sanity check REF
    int i, nref, reflen = strlen(line->d.allele[0]);
    char *ref = faidx_fetch_seq(args->fai, (char *) hdr->id[BCF_DT_CTG][line->rid].key, line->pos,
                                line->pos + reflen - 1, &nref);
    if (!ref)
    { error("faidx_fetch_seq failed at %s:%d\n", hdr->id[BCF_DT_CTG][line->rid].key, line->pos + 1); }
    replace_iupac_codes(ref, nref);


    // does REF contain non-standard bases?
    if (replace_iupac_codes(line->d.allele[0], reflen))
    {
        args->nchanged++;
        bcf_update_alleles(hdr, line, (const char **) line->d.allele, line->n_allele);
    }
    if (strcasecmp(ref, line->d.allele[0]))
    {
        // we will handle erros within the Normaliser class - jared
//        if (args->check_ref == CHECK_REF_EXIT)
//        {
//            error("Reference allele mismatch at %s:%d .. REF_SEQ:'%s' vs VCF:'%s'\n", bcf_seqname(hdr, line),
//                  line->pos + 1, ref, line->d.allele[0]);
//        }
//        if (args->check_ref & CHECK_REF_WARN)
//        {
//            fprintf(stderr, "REF_MISMATCH\t%s\t%d\t%s\n", bcf_seqname(hdr, line), line->pos + 1,
//                    line->d.allele[0]);
//        }
        free(ref);
        return ERR_REF_MISMATCH;
    }
    free(ref);
    ref = NULL;

    if (line->n_allele == 1)
    { return ERR_OK; }    // a REF

    // make a copy of each allele for trimming
    hts_expand0(kstring_t, line->n_allele, args->ntmp_als, args->tmp_als);
    kstring_t *als = args->tmp_als;
    for (i = 0; i < line->n_allele; i++)
    {
        if (line->d.allele[i][0] == '<') return ERR_SYMBOLIC;  // symbolic allele

        als[i].l = 0;
        kputs(line->d.allele[i], &als[i]);

        if (i > 0 && als[i].l == als[0].l && !strcasecmp(als[0].s, als[i].s)) return ERR_DUP_ALLELE;
    }


    // trim from right
    int ori_pos = line->pos;
    while (1)
    {
        // is the rightmost base identical in all alleles?
        for (i = 1; i < line->n_allele; i++)
        {
            if (als[0].s[als[0].l - 1] != als[i].s[als[i].l - 1]) break;
        }
        if (i != line->n_allele) break; // there are differences, cannot be trimmed

        int pad_from_left = 0;
        for (i = 0; i < line->n_allele; i++) // trim all alleles
        {
            als[i].l--;
            if (!als[i].l) pad_from_left = 1;
        }
        if (pad_from_left)
        {
            int npad = line->pos >= args->aln_win ? args->aln_win : line->pos;
            free(ref);
            ref = faidx_fetch_seq(args->fai, (char *) hdr->id[BCF_DT_CTG][line->rid].key, line->pos - npad,
                                  line->pos - 1, &nref);
            if (!ref)
                error("faidx_fetch_seq failed at %s:%d\n", hdr->id[BCF_DT_CTG][line->rid].key,
                      line->pos - npad + 1);
            replace_iupac_codes(ref, nref);
            for (i = 0; i < line->n_allele; i++)
            {
                ks_resize(&als[i], als[i].l + npad);
                if (als[i].l) memmove(als[i].s + npad, als[i].s, als[i].l);
                memcpy(als[i].s, ref, npad);
                als[i].l += npad;
            }
            line->pos -= npad;
        }
    }
    free(ref);

    // trim from left
    int ntrim_left = 0;
    while (1)
    {
        // is the first base identical in all alleles?
        int min_len = als[0].l - ntrim_left;
        for (i = 1; i < line->n_allele; i++)
        {
            if (als[0].s[ntrim_left] != als[i].s[ntrim_left]) break;
            if (min_len > als[i].l - ntrim_left) min_len = als[i].l - ntrim_left;
        }
        if (i != line->n_allele || min_len == 1) break; // there are differences, cannot be trimmed
        ntrim_left++;
    }
    if (ntrim_left)
    {
        for (i = 0; i < line->n_allele; i++)
        {
            memmove(als[i].s, als[i].s + ntrim_left, als[i].l - ntrim_left);
            als[i].l -= ntrim_left;
        }
        line->pos += ntrim_left;
    }

    // Have the alleles changed?
    als[0].s[als[0].l] = 0;  // in order for strcmp to work
    if (ori_pos == line->pos && !strcasecmp(line->d.allele[0], als[0].s)) return ERR_OK;

    // Create new block of alleles and update
    args->tmp_als_str.l = 0;
    for (i = 0; i < line->n_allele; i++)
    {
        if (i > 0) kputc(',', &args->tmp_als_str);
        kputsn(als[i].s, als[i].l, &args->tmp_als_str);
    }
    args->tmp_als_str.s[args->tmp_als_str.l] = 0;
    bcf_update_alleles_str(hdr, line, args->tmp_als_str.s);
    args->nchanged++;

    return ERR_OK;
}

static int diploid_to_haploid(int size, int nsmpl, int nals, uint8_t *vals)
{
    int i, dsrc = size * nals * (nals + 1) / 2, ddst = size * nals;
    uint8_t *src_ptr = vals + dsrc, *dst_ptr = vals + ddst;
    for (i = 1; i < nsmpl; i++)
    {
        memmove(dst_ptr, src_ptr, ddst);
        dst_ptr += ddst;
        src_ptr += dsrc;
    }
    return nals;
}

static void init_data(args_t *args)
{
    int nsample=1;
    args->hdr = NULL;
    rbuf_init(&args->rbuf, 100);
    args->lines = (bcf1_t **) calloc(args->rbuf.m, sizeof(bcf1_t *));
    if (args->ref_fname)
    {
        args->fai = fai_load(args->ref_fname);
        if (!args->fai)
        { error("Failed to load the fai index: %s\n", args->ref_fname); }
    }
    if (args->mrows_op == MROWS_MERGE)
    {
        args->mrow_out = bcf_init1();
        args->tmp_str = (kstring_t *) calloc(nsample, sizeof(kstring_t));
        args->diploid = (uint8_t *) malloc(nsample);
    }
    args->check_ref |= CHECK_REF_WARN;
}

void destroy_data(args_t *args)
{
    int i;
    for (i = 0; i < args->rbuf.m; i++)
        if (args->lines[i]) bcf_destroy1(args->lines[i]);
    free(args->lines);
    for (i = 0; i < args->mtmp_lines; i++)
        if (args->tmp_lines[i]) bcf_destroy1(args->tmp_lines[i]);
    free(args->tmp_lines);
    for (i = 0; i < args->malines; i++)
        bcf_destroy1(args->alines[i]);
    free(args->alines);
    for (i = 0; i < args->mblines; i++)
        bcf_destroy1(args->blines[i]);
    free(args->blines);
    for (i = 0; i < args->mmaps; i++)
        free(args->maps[i].map);
    for (i = 0; i < args->ntmp_als; i++)
        free(args->tmp_als[i].s);
    free(args->tmp_als);
    free(args->tmp_als_str.s);
    if (args->tmp_str)
    {
        for (i = 0; i < bcf_hdr_nsamples(args->hdr); i++) free(args->tmp_str[i].s);
        free(args->tmp_str);
    }
    free(args->maps);
    free(args->als);
    free(args->tmp_arr1);
    free(args->tmp_arr2);
    free(args->diploid);
    if (args->mrow_out) bcf_destroy1(args->mrow_out);
    if (args->fai)
    { fai_destroy(args->fai); }
    if (args->mseq)
    { free(args->seq); }
}

args_t *init_vcfnorm(bcf_hdr_t *hdr,  char *ref)
{
    args_t *args = (args_t *) calloc(1, sizeof(args_t));
    args->files = NULL;
    args->output_fname = "-";
    args->output_type = FT_VCF;
    args->aln_win = 100;
    args->buf_win = 1000;
    args->mrows_collapse = COLLAPSE_BOTH;
    args->mrows_op = MROWS_SPLIT;
    args->hdr = hdr;
    args->do_indels = 1;
    args->ref_fname = ref;
    init_data(args);
    return (args);
}
