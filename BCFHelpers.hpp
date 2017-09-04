// -*- mode: c++; indent-tabs-mode: nil; -*-
//
//
// Copyright (c) 2010-2015 Illumina, Inc.
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice,
// this
//    list of conditions and the following disclaimer.

// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

/**
 *  \brief C++ Helper functions for HTSlib
 *
 *
 * \file BCFHelpers.hpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include <iostream>
#include <list>
#include <memory>
#include <queue>
#include <vector>

extern "C" {

// GCC warns us about some things in htslib here. We don't care.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wredundant-decls"

#include <htslib/synced_bcf_reader.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>

#pragma GCC diagnostic pop
}

#include "Fasta.hpp"
// #include "Error.hpp"
#include "RefVar.hpp"

namespace bcfhelpers
{

extern const int MAX_GT;

/** when things go wrong, this gets thrown */
struct importexception : public std::runtime_error
{
    importexception(std::string const what)
            : std::runtime_error(what)
    {
    }
};

/** make a missing float */
static float missing_float() __attribute__((unused));

static float missing_float()
{
    union
    {
        int32_t i;
        float f;
    } val;
    val.i = bcf_float_missing;
    return val.f;
}

/* get list of sample names from BCF header */
std::list<std::string> getSampleNames(const bcf_hdr_t *hdr);

/** shortcut to get chromosome name */
static inline std::string getChrom(const bcf_hdr_t *hdr, const bcf1_t *rec)
{
    return hdr->id[BCF_DT_CTG][rec->rid].key;
}

/** extract pos / length for all samples in record */
void getLocation(bcf_hdr_t *hdr, bcf1_t *rec, int64_t &start, int64_t &end);

/**
 * @brief Retrieve an info field as an integer
 * @details [long description]
 *
 * @param result the default to return if the field is not present
 */
std::string getInfoString(bcf_hdr_t *header, bcf1_t *line, const char *field, const char *def_result = ".");

/**
 * @brief Retrieve an info field as an integer
 * @details [long description]
 *
 * @param result the default to return if the field is not present
 */
int getInfoInt(bcf_hdr_t *header, bcf1_t *line, const char *field, int result = -1);

std::vector<int> getInfoInts(bcf_hdr_t *header, bcf1_t *line, const char *field);

/**
 * @brief Retrieve an info field as a double
 *
 * @return the value or NaN
 */
float getInfoFloat(bcf_hdr_t *header, bcf1_t *line, const char *field);

std::vector<float> getInfoFloats(bcf_hdr_t *header, bcf1_t *line, const char *field);

/**
 * @brief Retrieve an info flag
 *
 * @return true of the flag is set
 */
bool getInfoFlag(bcf_hdr_t *header, bcf1_t *line, const char *field);

/**
 * @brief Read the GT field
 */
void getGT(bcf_hdr_t *header, bcf1_t *line, int isample, int *gt, int &ngt, bool &phased);

/** read GQ(X) -- will use in this order: GQ, GQX, -1 */
void getGQ(const bcf_hdr_t *header, bcf1_t *line, int isample, float &gq);

/** read AD */
void getAD(const bcf_hdr_t *header, bcf1_t *line, int isample, int *ad, int max_ad);

/** read DP(I) -- will use in this order: DP, DPI, -1 */
void getDP(const bcf_hdr_t *header, bcf1_t *line, int isample, int &dp);

/** read a format field as a single int. */
int getFormatInt(const bcf_hdr_t *header, bcf1_t *line, const char *field, int isample, int defaultresult = -1);

std::vector<int> getFormatInts(const bcf_hdr_t *header, bcf1_t *line, const char *field, int isample);

/** read a format field as a single double. default return value is NaN */
float getFormatFloat(const bcf_hdr_t *header, bcf1_t *line, const char *field, int isample);

std::vector<float> getFormatFloats(const bcf_hdr_t *header, bcf1_t *line, const char *field, int isample);

/** read a format field as string. result will not be overwritten on failure */
std::string
getFormatString(const bcf_hdr_t *header, bcf1_t *line, const char *field, int isample, const char *def_result = ".");

/** update format string for a single sample.  */
void setFormatStrings(const bcf_hdr_t *header, bcf1_t *line, const char *field, const std::vector<std::string> &value);

/** update format with single float values.  */
void setFormatFloats(const bcf_hdr_t *header, bcf1_t *line, const char *field, const std::vector<float> &value);

/** update format with single int values.  */
void setFormatInts(const bcf_hdr_t *header, bcf1_t *line, const char *field, const std::vector<int> &value);

/** return number of reference padding bases */
int isRefPadded(bcf1_t *line);

/** return number of reference padding bases */
int addRefPad(bcf_hdr_t *hdr, bcf1_t *line, FastaFile const &ref, int npad = 1);

/** check if a variant is a (possibly multiallelic) SNP */
bool isSNP(bcf_hdr_t *hdr, bcf1_t *rec);

/**
 * Check if variant is a single-allelic MNP ie. can be trivially decomposed
 */
bool isMNP(bcf_hdr_t *hdr, bcf1_t *rec);

/**
 * Check if variant has any complex alleles which are neither insertions nor deletions
 *
 * This function detects reference padding.
 */
bool isComplex(bcf_hdr_t *hdr, bcf1_t *rec);

/**
 * Check if a variant has symbolic alleles beyond deletions and overlapping
 * or marginalized alleles
 */
bool isSymbolic(bcf_hdr_t *hdr, bcf1_t *rec);

/** Check if variant is an indel or complex event */
inline bool isINDEL(bcf_hdr_t *hdr, bcf1_t *rec)
{
    bcf_unpack(rec, BCF_UN_SHR);
    if(rec->n_allele == 1)
    {
        return false;
    }
    if(isSNP(hdr, rec))
    {
        return false;
    }
    if(isMNP(hdr, rec))
    {
        return false;
    }
    if(isSymbolic(hdr, rec))
    {
        return false;
    }
    return (true);
}

/**
 * @brief adds an ALT=<M> variant to bcf record
 *
 * @param rec a points to a bcf1_t
 * @param hdr the respective header for rec
 *
 */
void addSymbolic(const bcf_hdr_t *hdr, bcf1_t *rec);

/** shared pointer support for keeping bcf types around */
typedef std::shared_ptr<bcf_srs_t> p_bcf_srs_t;
typedef std::shared_ptr<bcf_hdr_t> p_bcf_hdr;
typedef std::shared_ptr<bcf1_t> p_bcf1;
typedef std::shared_ptr<htsFile> p_hts_file;
typedef std::shared_ptr<hts_idx_t> p_hts_index;
typedef std::shared_ptr<tbx_t> p_tbx;

/** shared-pointerize */
static inline p_bcf_hdr p(bcf_hdr_t *h)
{
    return std::shared_ptr<bcf_hdr_t>(h, bcf_hdr_destroy);
}

static inline p_bcf1 p(bcf1_t *b)
{
    return std::shared_ptr<bcf1_t>(b, bcf_destroy);
}

static inline p_bcf_srs_t p(bcf_srs_t *b)
{
    return std::shared_ptr<bcf_srs_t>(b, bcf_sr_destroy);
}

static inline p_hts_file p(htsFile *b)
{
    return std::shared_ptr<htsFile>(b, hts_close);
}

static inline p_hts_index p(hts_idx_t *b)
{
    return std::shared_ptr<hts_idx_t>(b, hts_idx_destroy);
}

static inline p_tbx p(tbx_t *b)
{
    return std::shared_ptr<tbx_t>(b, tbx_destroy);
}

/** split multi-allelic variants in line into bi-allelic variants
 * this is modifed version of bcftools norm function
 *
 * We return data in a vector of shared pointers, this is to make clear
 * who owns the data.
 *
 */
void splitMultiAllelics(bcf_hdr_t *hdr, bcf1_t *line, std::vector<p_bcf1> &out);

/**
 * Compare bcf records by position
 * @param hdr bcf header
 * @param a first bcf record
 * @param b second bcf record
 * @return true if a < b else false
 */
bool inline compare(bcf_hdr_t *hdr, bcf1_t *a, bcf1_t *b)
{
    if(a->pos == b->pos)
    {
        int *gt = NULL;
        int ngt = 0;
        int na = bcf_get_genotypes(hdr, a, &gt, &ngt);
        int nb = bcf_get_genotypes(hdr, b, &gt, &ngt);
        free(gt);
        return (a->n_allele < b->n_allele && na >= nb);
    }
    else
    {
        return (a->pos < b->pos);
    }
}

/**
 * Check bcf records for equality
 * Equality is determined via chrom, pos, alleles
 * @param a first bcf record
 * @param b second bcf record
 * @return true if a == b
 */
bool inline equal(bcf1_t *a, bcf1_t *b)
{
    if(a->rid != b->rid)
    {
        return false;
    }
    if(a->pos != b->pos || a->n_allele != b->n_allele)
    {
        return false;
    }
    for(int i = 0; i < a->n_allele; i++)
    {
        if(strcmp(a->d.allele[i], b->d.allele[i]) != 0)
        {
            return false;
        }
    }
    return true;
}

/**
 * @brief produces a string represent the variant (for printing etc)
 *
 * @param rec a points to a bcf1_t
 * @param hdr the respective header for rec
 *
 * @return a string "CHR:POS:REF/ALT"
 *
 */
static inline std::string summarize(bcf_hdr_t *hdr, bcf1_t *rec)
{
    bcf_unpack(rec, BCF_UN_SHR);
    std::string ret = bcfhelpers::getChrom(hdr, rec) + ":" + std::to_string(rec->pos + 1) + ":" + rec->d.allele[0];
    for(int i = 1; i < rec->n_allele; i++)
    {
        ret += "/" + (std::string) rec->d.allele[i];
    }
    return (ret);
}

static inline bool failFilter(bcf_hdr_t *hdr, bcf1_t *rec)
{
    bool fail = false;
    bcf_unpack(rec, BCF_UN_FLT);
    for(int j = 0; j < rec->d.n_flt; ++j)
    {
        std::string filter = "PASS";
        int k = rec->d.flt[j];

        if(k >= 0)
        {
            filter = bcf_hdr_int2id(hdr, BCF_DT_ID, rec->d.flt[j]);
        }
        if(filter != "PASS")
        {
            fail = true;
            break;
        }
    }
    return fail;
}

static inline bool failFilter(bcf_hdr_t *hdr, bcf1_t *rec,int sample_index)
{
    bcf_unpack(rec, BCF_UN_ALL);
    std::string sample_filter = bcfhelpers::getFormatString(hdr, rec, "FT", sample_index, ".");
    return sample_filter.compare("PASS");
}

typedef std::deque<p_bcf1> variantqueue_t;
enum class AlleleType
{
    NUC,              // Nucleotide allele
    MISSING,          // no allele / .
    SYMBOLIC_DEL,     // <DEL>
    OVERLAPPING_DEL,  // *
    MARGINALIZED,     // <M> / <*>
    SYMBOLIC_OTHER,   // Any other symbolic allele
    UNKNOWN           // unknown allele that is not symbolic
};

std::ostream &operator<<(std::ostream &, const AlleleType);

/**
 * Classify an allele string from VCF and transform / canonicalize some
 * representational differences ("" becomes ".", deletions become "", ...)
 * @param al the allele string as it comes from the BCF record
 * @return allele type and transformed / standardised / uppercased allele string
 */
std::pair<AlleleType, std::string> classifyAlleleString(std::string al);

/**
 * Produce AD record for a subset of the alleles, assuming that
 * the other alleles are matched to the reference (e.g. when complex
 * alleles get decomposed and produce a reference base).
 *
 * Requires:
 * - number of ALT alleles must be <= 64
 * - ref_alleles & keep_alleles = 0
 *
 * @param hdr BCF header
 * @param in BCF record
 * @param isample the sample to use
 * @param keep_alleles mask for the alt alleles
 * @param ref_alleles mask for alleles that are now REF
 * @param new_ad vector that will be resized to popcount(keep_alleles) + 1, gives the new ADs
 * @param marginalize final allele in new_ad will be <M> allele with depth over all non-ref, non-kept ALs
 */
void decomposeAD(bcf_hdr_t *hdr,
                 bcf1_t *in,
                 int isample,
                 uint64_t keep_alleles,
                 uint64_t ref_alleles,
                 std::vector<int> &new_ad,
                 bool marginalize = false);

/**
 * Create allele masks which specify for each PL which allele they correspond
 * to. The LSB will denote the REF allele (so to use these masks in any decompose
 * function, right-shift by one).
 *
 * This function follows the VCF Spec https://samtools.github.io/hts-specs/VCFv4.3.pdf
 * (page 10)
 *
 * @param ngt number of genotype entries (ploidy)
 * @param nal number of different alleles including reference
 * @return a vector of lists which give the GT alleles in sorted order
 */
std::vector<std::list<int> > const &makePLAlleleLists(int ngt, int nal);

/**
 * Produce PL record for a subset of the alleles, assuming that
 * the other alleles are matched to the reference (e.g. when complex
 * alleles get decomposed and produce a reference base).
 *
 * Requires:
 * - number of ALT alleles must be <= 64
 * - ref_alleles & keep_alleles = 0
 * - ploidy is 1 or 2
 *
 * @param hdr BCF header
 * @param in BCF record
 * @param isample the sample to use
 * @param keep_alleles mask for the alt alleles
 * @param ref_alleles mask for alleles that are now REF
 * @param new_pl updated PL values
 * @param marginalize final allele in new_ad will be <M> allele with marginalized PL for non-kept, non-ref alleles
 */
void decomposePL(bcf_hdr_t *hdr,
                 bcf1_t *in,
                 int isample,
                 uint64_t keep_alleles,
                 uint64_t ref_alleles,
                 std::vector<int> &new_pl,
                 bool marginalize = false);

/**
 * Retrieve GL values from a record, from PL if available
 *
 * Also checks the number of values is correct
 *
 * @param hdr BCF header
 * @param in record to read from
 * @param gl target array
 */
void getGL(bcf_hdr_t *hdr,
           bcf1_t *in,
           std::vector<float> &gl);

/**
 * Update CIGAR INFO field (if it exists)
 * @param hdr bcf header
 * @param rec bcf record
 */
void updateCigar(bcf_hdr_t *hdr, bcf1_t *rec);

} // namespace bcfhelpers
