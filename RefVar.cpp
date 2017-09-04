// -*- mode: c++; indent-tabs-mode: nil; -*-
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
 * \brief Reference Variation helper implementation
 *
 *
 * \file RefVar.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "RefVar.hh"

#include "Error.hh"
#include "Genetics.hh"

#include <algorithm>
#include <cassert>
#include <limits>

using namespace genetics;

// #define DEBUG_REFVAR

namespace variant
{

void trimLeft(FastaFile const &f, const char *chr, RefVar &rv, bool refpadding)
{
    // trim left
    std::string ref = f.query(chr, rv.start, rv.end);
    int64_t rel_start = 0;
    size_t ref_min = refpadding ? 1 : 0;

    while(ref.size() - rel_start > ref_min && rv.alt.size() - rel_start > ref_min &&
          ref[rel_start] == rv.alt[rel_start])
    {
        rel_start++;
        rv.start++;
    }
    if(rel_start > 0)
    {
        rv.alt = rv.alt.substr((unsigned long) rel_start);
    }
}

void trimRight(FastaFile const &f, const char *chr, RefVar &rv, bool refpadding)
{
    // trim right

    int64_t reflen = rv.end - rv.start + 1;
    int64_t altlen = (int64_t) rv.alt.size();
    int64_t min_len = refpadding ? 1 : 0;

    if(reflen <= min_len || altlen <= min_len)
    {
        return;
    }

    std::string ref = f.query(chr, rv.start, rv.end);

#ifdef DEBUG_REFVAR
    std::cerr << rv << " -- ref sq = " << ref << "\n";
#endif

    while(reflen > min_len && altlen > min_len && ref[reflen - 1] == rv.alt[altlen - 1])
    {
        altlen--;
        reflen--;
    }
    rv.end = rv.start + reflen - 1;
    if(altlen > 0)
    {
        rv.alt = rv.alt.substr(0, (unsigned long) altlen);
    }
    else
    {
        rv.alt = "";
    }
}

void leftShift(FastaFile const &f, const char *chr, RefVar &rv, int64_t pos_min)
{
    int64_t rstart = -1, rend = -1, reflen;

    // vaguely like
    // http://genome.sph.umich.edu/wiki/File:Variant_normalization_algorithm.png
    bool done = false;
    pos_min = std::max(pos_min, (int64_t) 0);

    trimLeft(f, chr, rv);
    trimRight(f, chr, rv);

    reflen = rv.end - rv.start + 1;
    // check for all ref match (HAP-64)
    if(reflen < 0 && rv.alt.size() == 0)
    {
        // no inserted allele and ref length < 0
        return;
    }

    std::string ref;
    if(reflen >= 0 && reflen == (signed) rv.alt.size())
    {
        ref = f.query(chr, rv.start, rv.end);
        if(ref == rv.alt)
        {
            return;
        }
    }

    while(!done)
    {
        done = true;
        reflen = rv.end - rv.start + 1;

        if(rstart < 0 || rv.start <= rstart || rend < 0 || rv.end > rend)
        {
            rstart = rv.start - 20;
            if(reflen <= 0)
            {
                rend = rv.start;
            }
            else
            {
                rend = rv.end;
            }

            if(rstart < 0)
            {
                rstart = 0;
            }
            ref = f.query(chr, rstart, rend);
        }
        if(rv.start <= pos_min)
        {
            break;
        }

        int64_t rel_start = rv.start - rstart;

        if(rel_start < 1 || ref.size() == 0 || ((signed) ref.size()) < rel_start + reflen ||
           ref[rel_start - 1] == 'N') // don't shift over Ns
        {
            break;
        }

        // right trim.
        if(reflen > 0 && rv.alt.size() > 0 && ref[rel_start + reflen - 1] == rv.alt.back())
        {
            reflen--;
            rv.end--;
            rv.alt = rv.alt.substr(0, rv.alt.size() - 1);
            done = false;
        }

        if(reflen == 0 || rv.alt.size() == 0)
        {
            rv.start--;
            rv.alt = ref.substr((unsigned long) (rel_start - 1), 1) + rv.alt;
            done = false;
        }
    }
    trimLeft(f, chr, rv);
    trimRight(f, chr, rv);
}

void rightShift(FastaFile const &f, const char *chr, RefVar &rv, int64_t pos_max)
{
    int64_t rstart = -1, rend = -1, reflen;

    trimLeft(f, chr, rv);
    trimRight(f, chr, rv);

    reflen = rv.end - rv.start + 1;
    // check for all ref match (HAP-64)
    if(reflen < 0 && rv.alt.size() == 0)
    {
        // no inserted allele and ref length < 0
        return;
    }

    if(reflen >= 0 && reflen == (signed) rv.alt.size())
    {
        std::string ref = f.query(chr, rv.start, rv.end);
        if(ref == rv.alt)
        {
            return;
        }
    }

    // adapted from
    // http://genome.sph.umich.edu/wiki/File:Variant_normalization_algorithm.png
    bool done = false;
    std::string ref;
    while(!done)
    {
        done = true;
        reflen = rv.end - rv.start + 1;

        if(rstart < 0 || rv.start < rstart || rend < 0 || rv.end >= rend)
        {
            rstart = rv.start;
            if(reflen <= 0)
            {
                rend = rv.start + 20;
            }
            else
            {
                rend = rv.end + 20;
            }

            if(rstart < 0)
            {
                rstart = 0;
            }
            ref = f.query(chr, rstart, rend);
        }
        if(rv.end >= pos_max)
        {
            break;
        }

        int64_t rel_start = rv.start - rstart;

        if(ref.size() == 0 || ((signed) ref.size()) <= rel_start + reflen || ref[rel_start + reflen] == 'N')
        {
            break;
        }

        // left trim.
        if(reflen > 0 && rv.alt.size() > 0 && ref[rel_start] == rv.alt[0])
        {
            reflen--;
            rel_start++;
            rv.start++;
            rv.alt = rv.alt.substr(1, rv.alt.size());
            done = false;
        }

        // right-extend
        if(reflen == 0 || rv.alt.size() == 0)
        {
            int64_t refnext = rel_start + reflen;
            rv.end++;
            rv.alt = rv.alt + ref.substr((unsigned long) refnext, 1);
            done = false;
        }
    }

    trimLeft(f, chr, rv);
    trimRight(f, chr, rv);
}

/**
 * Convert a list of RefVar records to allele strings
 */
void toAlleles(FastaFile const &f, const char *chr, std::vector<RefVar> const &in, std::vector<std::string> &out)
{
    int64_t minpos = std::numeric_limits<int64_t>::max(), maxpos = 0;

    if(in.size() == 0)
    {
        return;
    }

    for(size_t s = 0; s < in.size(); ++s)
    {
        RefVar const &rv = in[s];
        minpos = std::min(minpos, rv.start);
        minpos = std::min(minpos, rv.end);
        // insertions
        if(rv.start > rv.end)
        {
            maxpos = std::max(maxpos, rv.end);
        }
        else
        {
            maxpos = std::max(maxpos, rv.start);
            maxpos = std::max(maxpos, rv.end);
        }
    }
    std::string ref = f.query(chr, minpos, maxpos);
    if((signed) ref.size() != maxpos - minpos + 1)
    {
        error("Cannot query reference sequence at %s:%i-%i", ref.c_str(), minpos, maxpos);
    }
    out.resize(in.size() + 1);
    out[0] = ref;
    for(size_t i = 0; i < in.size(); ++i)
    {
        // ref allele for this one:
        // ref.substr(in[i].start - minpos, in[i].end - in[i].start + 1);
        int64_t refstart = in[i].start - minpos;
        int64_t reflen = in[i].end - in[i].start + 1;
        // don't break symbolic alleles
        if(bcfhelpers::classifyAlleleString(in[i].alt).first != bcfhelpers::AlleleType::NUC ||
           (reflen == (signed) ref.size() && refstart == 0))
        {
            out[i + 1] = in[i].alt;
        }
        else
        {
            // create alt allele by replacing
            out[i + 1] = ref;
            out[i + 1].replace((unsigned long) refstart, (unsigned long) reflen, in[i].alt);
        }
        if(out[i + 1].size() == 0)
        {
            out[i + 1] = "<DEL>";
        }
    }
}

/**
 * Make RefVar
 */
void mkRefVar(int64_t refpos, char refchr, char altchr, RefVar &rv)
{
    rv.start = refpos;
    rv.end = refpos;
    if(refchr == '-')
    {
        // alt insertion, refchr != altchr so altchr != '-'
        rv.end--;
        rv.alt = altchr;
    }
    else if(altchr == '-')
    {
        // ref deletion
        rv.alt = "";
    }
    else
    {
        // snp
        rv.alt = altchr;
    }
}

/**
 * Append 1 char to RefVar.
 *
 * Return false if this is not possible because the variants are incompatible
 */
bool appendToRefVar(int64_t refpos, char refchr, char altchr, RefVar &rv)
{
#ifdef DEBUG_REFVAR
    std::cerr << "\t" << refpos << " " << refchr << " " << altchr << " " << rv
              << "\n";
#endif
    if(refchr == altchr)
    {
        // already checked in function just below actually
        // we test anyway in case we ever want to expose this
        // function
        return false;
    }
    // see what we have currently
    if((refpos == rv.start - 1 || refpos == rv.start) && rv.end == rv.start - 1 &&
       rv.alt.size() > 0) // rv is an insertion just before refpos, or a complex
        // insertion ending at refpos
    {
        // we've been passed an insertion and we can extend?
        if(refchr == '-') //  -> && altchr != '-' since refchr != altchr
        {
            rv.alt += altchr;
        }
            // we've been passed a substitution?
        else if(altchr != '-')
        {
            if(rv.start > refpos)
            {
                // first subst after insertion -- fix start and end
                rv.start = refpos;
                rv.end = refpos;
            }
            else
            {
                rv.end++;
            }
            rv.alt += altchr;
        }
        else   // -> altchr == '-' && refchr != '-'
            //  above otherwise
        {
            // we've been passed a deletion
            if(rv.start > refpos)
            {
                // first deletion after insertion -- fix start and end
                rv.start = refpos;
                rv.end = refpos;
            }
            else
            {
                rv.end++;
            }
        }
        return true;
    }

    // here, we append only
    if(refpos != rv.end + 1)
    {
        return false;
    }

    if(rv.end >= rv.start && rv.alt.size() > 0) // block substitution
    {
        // we've been passed a substitution?
        if(refchr != '-' && altchr != '-')
        {
            rv.end++;
            rv.alt += altchr;
        }
        else if(refchr == '-') //  -> && altchr != '-' since refchr != altchr
        {
            // we've been passed an insertion
            rv.alt += altchr;
        }
        else // (altchr == '-') //  -> && refchr != '-' since refchr != altchr
        {
            // we've been passed a deletion
            rv.end++;
        }
    }
    else if(rv.end >= rv.start && rv.alt.size() == 0) // deletion
    {
        // we've been passed a substitution?
        if(refchr != '-' && altchr != '-')
        {
            rv.end++;
            rv.alt += altchr;
        }
        else if(refchr == '-') //  -> && altchr != '-' since refchr != altchr
        {
            // we've been passed an insertion
            rv.alt += altchr;
        }
        else  // if (altchr == '-') //  -> && refchr != '-' since refchr != altchr
        {
            // we've been passed a deletion
            rv.end++;
        }
    }
    return true;
}

/**
 * @brief Decompose a RefVar into primitive variants (subst / ins / del)
 *
 * @param rv the RefVar record
 * @param vars the primitive records
 */
int toPrimitives(FastaFile const &f, const char *chr, RefVar const &rv, std::list<variant::RefVar> &vars)
{
    int64_t rstart = rv.start, rend = rv.end, reflen = rend - rstart + 1;
    int64_t altlen = (int64_t) rv.alt.size();

    std::string refseq;
    std::string altseq(rv.alt);

    if(reflen > 0)
    {
        refseq = f.query(chr, rstart, rend);
    }

    // from the left, split off SNPs / matches
    size_t pos = 0;
    int count = 0;
    while(reflen > 0 && altlen > 0)
    {
        char r = refseq[pos];
        char a = altseq[pos];
        if(r != a)
        {
            vars.push_back(RefVar{rstart, rstart, std::string(1, a), rv.flags});
            count++;
        }
        ++rstart;
        ++pos;
        --reflen;
        --altlen;
    }
    // now either reflen == 0 or altlen == 0
    if(reflen > 0)
    {
        // del
        vars.push_back(RefVar{rstart, rend, "", rv.flags});
        count++;
    }
    else if(altlen > 0)
    {
        // ins
        vars.push_back(RefVar{rstart, rstart - 1, altseq.substr(pos), rv.flags});
        count++;
    }
    // else nothing left
    return (count);
}

/**
 * @brief Decompose a RefVar into primitive variants (subst / ins / del) by
 * means of realigning
 *
 * @param f reference sequence fasta
 * @param chr the chromosome to use
 * @param rv the RefVar record
 * @param snps the number of snps
 * @param ins the number of insertions
 * @param dels the number of deletions
 * @param homref the number of calls with no variation
 */
void countRefVarPrimitives(FastaFile const &f, const char *chr, variant::RefVar const &rv, size_t &snps, size_t &ins,
                           size_t &dels, size_t &homref, size_t &transitions, size_t &transversions)
{
    int64_t rstart = rv.start, rend = rv.end, reflen = rend - rstart + 1;
    int64_t altlen = (int64_t) rv.alt.size();

    std::string refseq;
    std::string altseq(rv.alt);

    if(reflen <= 0)
    {
        if(altlen > 0)
        {
            ins += altlen;
        }
        else
        {
            ++homref;
        }
        return;
    }
    // reflen > 0
    refseq = f.query(chr, rstart, rend);

    // from the left, split off SNPs / matches
    size_t pos = 0;
    bool isValidSnv(false);

    while(reflen > 0 && altlen > 0)
    {
        char r = refseq[pos];
        char a = altseq[pos];
        if(r != a)
        {
            ++snps;
            const bool isTransversion(snvIsTransversion(r, a, isValidSnv));

            if(isValidSnv)
            {
                if(isTransversion)
                {
                    ++transversions;
                }
                else
                {
                    ++transitions;
                }
            }
        }
        else
        {
            ++homref;
        }
        ++rstart;
        ++pos;
        --reflen;
        --altlen;
    }
    // now either reflen == 0 or altlen == 0
    if(reflen > 0)
    {
        // del
        dels += reflen;
    }
    if(altlen > 0)
    {
        // ins
        ins += altlen;
    }
    // else nothing left
}

/**
 * Helper to aggregate reference variants base by base
 *
 * Refpos must be > vars.back().end
 *
 */
void appendToVarList(int64_t refpos, char refchr, char altchr, std::list<variant::RefVar> &vars, int flags)
{
#ifdef DEBUG_REFVAR
    std::cerr << refpos << " " << refchr << " " << altchr
              << " nv: " << vars.size() << "\n";
#endif
    if(refchr == altchr)
    {
        return;
    }
    if(vars.empty())
    {
        RefVar rv;
        mkRefVar(refpos, refchr, altchr, rv);
        rv.flags = flags;
        vars.push_back(rv);
    }
    else
    {
        RefVar &xrv = vars.back();
        if(flags != xrv.flags || !appendToRefVar(refpos, refchr, altchr, xrv))
        {
            RefVar rv;
            mkRefVar(refpos, refchr, altchr, rv);
            rv.flags = flags;
            vars.push_back(rv);
        }
    }
}

/**
 * Trim multiple RefVars together
 * @param f a fasta reference
 * @param chr chromosme name
 * @param rv list of refvars
 * @param refpadding include reference padding
 */
int multiTrimLeft(FastaFile const &f, const char *chr, std::vector<variant::RefVar> &rv, bool refpadding)
{
    // make sure all variants are over the same interval.
    std::vector<bool> used(rv.size());
    for(size_t i = 0; i < rv.size(); i++)
    {
        auto type = bcfhelpers::classifyAlleleString(rv[i].alt).first;
        if(type == bcfhelpers::AlleleType::NUC)
        {
            used[i] = true;
        }
        else if(type == bcfhelpers::AlleleType::OVERLAPPING_DEL || type == bcfhelpers::AlleleType::MARGINALIZED)
        {
            used[i] = false;
        }
        else
        {
            return 0;
        }

        if(i > 0)
        {
            // make sure all variants are over the same interval.
            assert(rv[i].start == rv[i - 1].start && rv[i].end == rv[i - 1].end);
        }
    }
    int old_start = (int) rv[0].start;

    // trim left
    std::string ref = f.query(chr, rv[0].start, rv[0].end);
    int64_t rel_start = 0;
    size_t ref_min = refpadding ? 1 : 0;

    size_t min_alt_size = (size_t) (getMinAltLen(rv) - rel_start);
    bool start_equal = (ref[rel_start] == rv[0].alt[rel_start]);

    for(size_t i = 1; i < rv.size(); i++)
    {
        if(used[i])
        {
            start_equal = start_equal && (ref[rel_start] == rv[i].alt[rel_start]);
        }
    }

    while(ref.size() - rel_start > ref_min && min_alt_size > ref_min && start_equal)
    {
        rel_start++;
        for(size_t i = 0; i < rv.size(); i++)
        {
            rv[i].start++;
        }
        min_alt_size = (size_t) (getMinAltLen(rv) - rel_start);
        start_equal = (ref[rel_start] == rv[0].alt[rel_start]);
        for(size_t i = 1; i < rv.size(); i++)
        {
            if(used[i])
            {
                start_equal = start_equal && (ref[rel_start] == rv[i].alt[rel_start]);
            }
        }
    }

    if(rel_start > 0)
    {
        for(size_t i = 0; i < rv.size(); i++)
        {
            if(used[i])
            {
                rv[i].alt = rv[i].alt.substr((unsigned long) rel_start);
            }
        }
    }

    return (int) (rv[0].start - old_start);
}

/**
 * Trim multiple refvars on the right
 * @param f reference fasta file
 * @param chr chromosome name
 * @param rv list of refvars
 * @param refpadding include reference padding
 */
int multiTrimRight(FastaFile const &f, const char *chr, std::vector<variant::RefVar> &rv, bool refpadding)
{
    assert(!rv.empty());

    std::vector<bool> used(rv.size());
    for(size_t i = 0; i < rv.size(); i++)
    {
        auto type = bcfhelpers::classifyAlleleString(rv[i].alt).first;
        if(type == bcfhelpers::AlleleType::NUC)
        {
            used[i] = true;
        }
        else if(type == bcfhelpers::AlleleType::OVERLAPPING_DEL || type == bcfhelpers::AlleleType::MARGINALIZED)
        {
            used[i] = false;
        }
        else
        {
            return 0;
        }
        if(i > 0)
        {
            // make sure all variants are over the same interval.
            assert(rv[i].start == rv[i - 1].start && rv[i].end == rv[i - 1].end);
        }
    }

    int old_end = (int) rv[0].end;

    // trim right
    int64_t min_len = refpadding ? 1 : 0;
    int64_t reflen = rv[0].end - rv[0].start + 1;
    int64_t min_altlen = getMinAltLen(rv);

    if(reflen <= min_len || min_altlen <= min_len)
    {
        return 0;
    }

    // fetch ref + one base before
    std::string ref = f.query(chr, rv[0].start - 1, rv[0].end);
    if(ref.size() != (unsigned) (rv[0].end - rv[0].start + 2))
    {
        error("problem extracting reference sequence at %i", rv[0].start);
    }

#ifdef DEBUG_REFVAR
    std::cerr << rv << " -- ref sq = " << ref << "\n";
#endif

    bool end_equal = true;
    for(size_t i = 0; i < rv.size(); i++)
    {
        if(!used[i])
        {
            continue;
        }
        if(rv[i].alt.empty())
        {
            end_equal = false;
            break;
        }
        end_equal = end_equal && (ref.back() == rv[i].alt.back());
    }

    int ntrim = 0;
    while(reflen > min_len && min_altlen > min_len && end_equal)
    {
        ntrim++;
        reflen--;
        min_altlen = getMinAltLen(rv);
        end_equal = reflen < (signed) ref.size();
        for(size_t i = 0; i < rv.size(); i++)
        {
            if(!used[i])
            {
                continue;
            }
            const uint64_t altpos = rv[i].alt.size() - 1 - ntrim;
            if(altpos < rv[i].alt.size()) // if we have completely trimmed away an allele stop.
            {
                end_equal = false;
                break;
            }
            end_equal = end_equal && (ref[reflen] == rv[i].alt[altpos]);
        }
    }

    for(size_t i = 0; i < rv.size(); i++)
    {
        rv[i].end = rv[i].start + reflen - 1;

        if(!used[i])
        {
            continue;
        }

        if(min_altlen - ntrim > 0)
        {
            rv[i].alt = rv[i].alt.substr(0, rv[i].alt.size() - ntrim);
        }
        else
        {
            rv[i].alt = "";
        }
    }
    return (int) (old_end - rv[0].end);
}

/**
 * Shift multiple RefVars left together (for variant normalisation where
 * we want to preserve ref overlap)
 * @param f reference fasta file
 * @param chr reference contig
 * @param rv list of refvars
 * @param pos_min minimum position for shifting
 * @return bit index abc. a=left-shifted. b=left-trimmed. c=right-trimmed.
 */
std::vector<bool> multiLeftShift(FastaFile const &f, const char *chr, std::vector<variant::RefVar> &rv, int64_t pos_min)
{
    // bit index abc. a=left-shifted. b=left-trimmed. c=right-trimmed.
    std::vector<bool> return_value(3, 0);

    // checks all alleles have consistent REF
    for(size_t i = 1; i < rv.size(); i++)
    {
        assert(rv[i].start == rv[i - 1].start && rv[i].end == rv[i - 1].end);
    }

    int64_t rstart = -1, rend = -1, reflen;

    // vaguely like
    // http://genome.sph.umich.edu/wiki/File:Variant_normalization_algorithm.png
    bool done = false;
    std::string ref;

    pos_min = std::max(pos_min, (int64_t) pos_min);
    return_value[1] = multiTrimLeft(f, chr, rv) > 0;
    return_value[2] = multiTrimRight(f, chr, rv) > 0;

    reflen = rv[0].end - rv[0].start + 1;

    // check for all ref match (HAP-64)
    if(reflen < 0 && getMinAltLen(rv) == 0)
    {
        // no inserted allele and ref length < 0
        return return_value;
    }

    // this is checking if the ALT are all the same as the REF
    if(reflen > 0)
    {
        ref = f.query(chr, rv[0].start, rv[0].end);
        bool all_equal = true;
        for(size_t i = 0; i < rv.size(); i++)
        {
            all_equal = all_equal && (ref == rv[i].alt);
        }
        if(all_equal)
        {
            return return_value;
        }
    }

    std::vector<bool> used(rv.size());
    for(size_t i = 0; i < rv.size(); i++)
    {
        auto type = bcfhelpers::classifyAlleleString(rv[i].alt).first;
        if(type == bcfhelpers::AlleleType::NUC || type == bcfhelpers::AlleleType::SYMBOLIC_DEL)
        {
            used[i] = true;
        }
        else if(type == bcfhelpers::AlleleType::OVERLAPPING_DEL || type == bcfhelpers::AlleleType::MARGINALIZED)
        {
            used[i] = false;
        }
        else
        {
            return return_value;
        }
        if(i > 0)
        {
            // make sure all variants are over the same interval.
            assert(rv[i].start == rv[i - 1].start && rv[i].end == rv[i - 1].end);
        }
    }

    int old_start = (int) rv[0].start;
    while(!done)
    {
        done = true;
        int64_t ref_start = rv[0].start;
        int64_t ref_end = rv[0].end;
        reflen = ref_end - ref_start + 1;

        // check if leftshift limit
        if(ref_start <= pos_min)
        {
            break;
        }

        // retrieve new reference sequence
        if(rstart < 0 || ref_start <= rstart || rend < 0 || ref_end > rend)
        {
            rstart = ref_start - 20;
            if(reflen <= 0)
            {
                rend = ref_start;
            }
            else
            {
                rend = ref_end;
            }
            if(rstart < 0)
            {
                rstart = 0;
            }
            ref = f.query(chr, rstart, rend);
        }

        int64_t rel_start = ref_start - rstart;

        if(rel_start < 1 || ref.size() == 0 || ((signed) ref.size()) < rel_start + reflen ||
           ref[rel_start - 1] == 'N') // don't shift over Ns
        {
            break;
        }

        // right trim.
        if(reflen > 0 && getMinAltLen(rv) > 0)
        {
            bool all_equal = true;
            for(size_t i = 0; i < rv.size(); i++)
            {
                if(!used[i])
                {
                    continue;
                }
                all_equal = all_equal && (ref[rel_start + reflen - 1] == rv[i].alt.back());
            }

            if(all_equal)
            {
                reflen--;
                for(size_t i = 0; i < rv.size(); i++)
                {
                    rv[i].end--;

                    if(!used[i])
                    {
                        continue;
                    }
                    rv[i].alt = rv[i].alt.substr(0, rv[i].alt.size() - 1);
                }
                done = false;
            }
        }

        if(reflen == 0 || getMinAltLen(rv) == 0)
        {
            for(size_t i = 0; i < rv.size(); i++)
            {
                rv[i].start--;

                if(!used[i])
                {
                    continue;
                }
                rv[i].alt = ref.substr((unsigned long) (rel_start - 1), 1) + rv[i].alt;
            }
            done = false;
        }
    }

    return_value[0] = old_start > rv[0].start ? 1 : 0;
    return_value[1] = multiTrimLeft(f, chr, rv) > 0 || return_value[1] ? 1 : 0;
    return_value[2] = multiTrimRight(f, chr, rv) > 0 || return_value[2] ? 1 : 0;

    return return_value;
}

/**
 * Get the length of the shortest ALT allele for a vector of refVars
 * @param rv
 * @return
 */
int64_t getMinAltLen(std::vector<variant::RefVar> &rv)
{
    if(rv.empty())
    {
        return 0;
    }
    size_t min_altlen = (size_t) -1;
    for(size_t i = 0; i < rv.size(); i++)
    {
        auto type = bcfhelpers::classifyAlleleString(rv[i].alt).first;
        if(type == bcfhelpers::AlleleType::NUC)
        {
            if(min_altlen == ((size_t) -1) || rv[i].alt.size() < min_altlen)
            {
                min_altlen = rv[i].alt.size();
            }
        }
        else if(type == bcfhelpers::AlleleType::SYMBOLIC_DEL)
        {
            return 0;
        }
    }

    if(min_altlen == (size_t) -1)
    {
        return 0;
    }

    return ((int64_t) min_altlen);
}

} // namespace variant
