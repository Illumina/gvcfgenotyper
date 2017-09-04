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
 *  \brief Implementation of BCF file format helpers
 *
 * \file BCFHelpers.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "BCFHelpers.hh"
#include "StringUtil.hh"

#include "PopCount.hh"
#include "Error.hh"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <set>
#include <sstream>

#include <htslib/vcf.h>
#include <memory>
#include <unordered_map>
#include <Alignment.hh>

/**
 * @brief Helper to get out GT fields
 */

namespace bcfhelpers
{
const int MAX_GT = 5;

namespace _impl
{
/** C++-ified version of bcf_get_format_values */
template<typename target_type_t>
struct bcf_get_numeric_format
{
    bcf_get_numeric_format()
    {}

    /* return true when successful
*/
    void
    operator()(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, int isample, std::vector<target_type_t> &dest) const
    {
        dest.clear();
        int i;
        int tag_id = bcf_hdr_id2int(hdr, BCF_DT_ID, tag);

        if(!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, tag_id))
        {
            return;
        }

        int nsmpl = bcf_hdr_nsamples(hdr);
        if(isample >= nsmpl)
        {
            return;
        }

        bcf_unpack(line, BCF_UN_FMT);

        for(i = 0; i < line->n_fmt; i++)
        {
            if(line->d.fmt[i].id == tag_id)
            {
                break;
            }
        }

        if(i == line->n_fmt)
        {
            return;
        }

        bcf_fmt_t *fmt = &line->d.fmt[i];
        int type = fmt->type;
        dest.resize((unsigned long) fmt->n);

        for(i = 0; i < fmt->n; ++i)
        {
            if(type == BCF_BT_FLOAT)
            {
                static const auto make_missing_float = []() -> float
                {
                    float f;
                    int32_t _mfloat = 0x7F800001;
                    memcpy(&f, &_mfloat, 4);
                    return f;
                };
                static const float bcf_missing_float = make_missing_float();
                // cppcheck-suppress invalidPointerCast
                float res = ((float *) (fmt->p + isample * fmt->size))[i];
                if(res != bcf_missing_float)
                {
                    dest[i] = target_type_t(res);
                }
            }
            else if(type == BCF_BT_INT8)
            {
                int8_t r = ((int8_t *) (fmt->p + isample * fmt->size))[i];
                if(r == bcf_int8_vector_end)
                {
                    break;
                }
                if(r != bcf_int8_missing)
                {
                    dest[i] = target_type_t(r);
                }
            }
            else if(type == BCF_BT_INT16)
            {
                int16_t r = ((int16_t *) (fmt->p + isample * fmt->size))[i];
                if(r == bcf_int16_vector_end)
                {
                    break;
                }
                if(r != bcf_int16_missing)
                {
                    dest[i] = target_type_t(r);
                }
            }
            else if(type == BCF_BT_INT32)
            {
                int32_t r = ((int32_t *) (fmt->p + isample * fmt->size))[i];
                if(r == bcf_int32_vector_end)
                {
                    break;
                }
                if(r != bcf_int32_missing)
                {
                    dest[i] = target_type_t(r);
                }
            }
            else
            {
                std::cerr << "[W] string format field ignored when looking for numeric "
                        "formats!"
                          << "\n";
                dest[i] = -1;
            }
        }
    }
};

template<typename type_t>
struct bcf_get_gts
{

    bcf_get_gts()
    {}

    /**
*@brief Extract GT from Format field
*
*/
    void operator()(bcf_fmt_t *gt, int i, int *igt, int & ngt, bool &phased) const
    {
        phased = false;
        type_t *p = (type_t *) (gt->p + i * gt->size);
        int ial;
        for(ial = 0; ial < MAX_GT; ial++)
        {
            igt[ial] = -1;
        }
        for(ial = 0; ial < gt->n; ial++)
        {
            if(p[ial] == vector_end)
            {
                ngt = ial;
                break;
            } /* smaller ploidy */
            if(!(p[ial] >> 1) || p[ial] == missing)
            {
                continue;
            } /* missing allele */
            int al = (p[ial] >> 1) - 1;
            igt[ial] = al;
            phased = phased || ((p[ial] & 1) != 0);
        }
    }

    static const type_t missing;
    static const type_t vector_end;
};

template<>
const int8_t bcf_get_gts<int8_t>::missing(bcf_int8_missing);

template<>
const int8_t bcf_get_gts<int8_t>::vector_end(bcf_int8_vector_end);

template<>
const int16_t bcf_get_gts<int16_t>::missing(bcf_int16_missing);

template<>
const int16_t bcf_get_gts<int16_t>::vector_end(bcf_int16_vector_end);

template<>
const int32_t bcf_get_gts<int32_t>::missing(bcf_int32_missing);

template<>
const int32_t bcf_get_gts<int32_t>::vector_end(bcf_int32_vector_end);

template<typename type_t>
struct bcf_get_info
{
    bcf_get_info()
    {}

    type_t operator()(bcf_info_t *field, int which = 0) const;
};

template<>
int bcf_get_info<int>::operator()(bcf_info_t *field, int which) const
{
    if(field->type != BCF_BT_CHAR && field->len <= which)
    {
        error("Cannot extract int from non-scalar INFO field (len = %i, requested: "
                      "%i).",
              field->len, which);
    }
    void *dataptr = nullptr;
    if(which == 0)
    {
        dataptr = &field->v1;
    }
    else
    {
        dataptr = field->vptr;
    }
    assert(dataptr);
    switch(field->type)
    {
        case BCF_BT_NULL:
            return 0;
        case BCF_BT_INT8:
            return (int) *(((int8_t *) dataptr) + which);
        case BCF_BT_INT16:
            return (int) *(((int16_t *) dataptr) + which);
        case BCF_BT_INT32:
            return (int) *(((int32_t *) dataptr) + which);
        case BCF_BT_FLOAT:
            return (int) *(((float *) dataptr) + which);
        case BCF_BT_CHAR:
            if(which > 0)
            {
                error("Cannot extract int %i from string INFO field", which);
            }
            return atoi((const char *) field->vptr);
        default:
            break;
    }
    return -1;
}

template<>
float bcf_get_info<float>::operator()(bcf_info_t *field, int which) const
{
    switch(field->type)
    {
        case BCF_BT_NULL:
            return std::numeric_limits<float>::quiet_NaN();
        case BCF_BT_INT8:
        case BCF_BT_INT16:
        case BCF_BT_INT32:
        {
            const bcf_get_info<int> gi;
            return (float) gi(field, which);
        }
        case BCF_BT_FLOAT:
            if(field->len <= which)
            {
                error("Cannot extract int from non-scalar INFO field (len = %i, "
                              "requested: %i).",
                      field->len, which);
            }
            return *(((float *) field->vptr) + which);
        case BCF_BT_CHAR:
            if(which > 0)
            {
                error("Cannot extract int %i from string INFO field", which);
            }
            return (float) atof((const char *) field->vptr);
        default:
            break;
    }
    return std::numeric_limits<float>::quiet_NaN();
}

template<>
std::string bcf_get_info<std::string>::operator()(bcf_info_t *field, int) const
{
    char num[256];
    switch(field->type)
    {
        case BCF_BT_NULL:
            return std::string("NULL");
        case BCF_BT_INT8:
        case BCF_BT_INT16:
        case BCF_BT_INT32:
            snprintf(num, 256, "%i", field->v1.i);
            return std::string(num);
        case BCF_BT_FLOAT:
            snprintf(num, 256, "%g", field->v1.f);
            return std::string(num);
        case BCF_BT_CHAR:
            if(field->vptr && field->len > 0)
            {
                return std::string((const char *) field->vptr, (unsigned long) field->len);
            }
        default:
            break;
    }
    return std::string();
}
}

/** extract chrom / pos / length */
void getLocation(bcf_hdr_t *hdr, bcf1_t *rec, int64_t &refstart, int64_t &refend)
{
    bcf_unpack(rec, BCF_UN_STR);
    refstart = rec->pos;
    refend = refstart;

    int endfield = getInfoInt(hdr, rec, "END", -1);

    if(endfield > 0)
    {
        // if there is an end field, don't validate the ref allele
        refend = endfield - 1;
    }
    else
    {
        refend = (int64_t) (refstart + strlen(rec->d.allele[0]) - 1);
        if(strchr(rec->d.allele[0], '.') || strchr(rec->d.allele[0], '-'))
        {
            // length might be inaccurate now
            refend = refstart;
            throw bcfhelpers::importexception(
                    std::string("[W] Unsupported REF allele with undefined length: ") + std::string(rec->d.allele[0]));
        }
    }
}

/**
 * @brief Retrieve an info field as an integer
 *
 * @param result the default to return if the field is not present
 */
std::string getInfoString(bcf_hdr_t *header, bcf1_t *line, const char *field, const char *def_result)
{
    std::string result = def_result;
    bcf_info_t *info_ptr = bcf_get_info(header, line, field);
    if(info_ptr)
    {
        static const bcfhelpers::_impl::bcf_get_info<std::string> i;
        result = i(info_ptr);
    }
    return result;
}

/**
 * @brief Retrieve an info field as an integer
 *
 * @param result the default to return if the field is not present
 */
int getInfoInt(bcf_hdr_t *header, bcf1_t *line, const char *field, int result)
{
    bcf_unpack(line, BCF_UN_INFO);
    bcf_info_t *info_ptr = bcf_get_info(header, line, field);
    if(info_ptr)
    {
        static const bcfhelpers::_impl::bcf_get_info<int> i;
        result = i(info_ptr);
    }
    return result;
}

std::vector<int> getInfoInts(bcf_hdr_t *header, bcf1_t *line, const char *field)
{
    std::vector<int> result;
    bcf_unpack(line, BCF_UN_INFO);
    bcf_info_t *info_ptr = bcf_get_info(header, line, field);
    if(info_ptr)
    {
        static const bcfhelpers::_impl::bcf_get_info<int> i;
        for(int j = 0; j < info_ptr->len; ++j)
        {
            result.push_back(i(info_ptr, j));
        }
    }
    return result;
}

/**
 * @brief Retrieve an info field as a double
 *
 * @return the value or NaN
 */
float getInfoFloat(bcf_hdr_t *header, bcf1_t *line, const char *field)
{
    float result = std::numeric_limits<float>::quiet_NaN();
    bcf_unpack(line, BCF_UN_INFO);
    bcf_info_t *info_ptr = bcf_get_info(header, line, field);
    if(info_ptr)
    {
        static const bcfhelpers::_impl::bcf_get_info<float> i;
        result = i(info_ptr);
    }
    return result;
}

std::vector<float> getInfoFloats(bcf_hdr_t *header, bcf1_t *line, const char *field)
{
    std::vector<float> result;
    bcf_unpack(line, BCF_UN_INFO);
    bcf_info_t *info_ptr = bcf_get_info(header, line, field);
    if(info_ptr)
    {
        static const bcfhelpers::_impl::bcf_get_info<float> i;
        for(int j = 0; j < info_ptr->len; ++j)
        {
            result.push_back(i(info_ptr, j));
        }
    }
    return result;
}

/**
 * @brief Retrieve an info flag
 *
 * @return true of the flag is set
 */
bool getInfoFlag(bcf_hdr_t *hdr, bcf1_t *line, const char *field)
{
    bcf_unpack(line, BCF_UN_INFO);
    return bcf_get_info_flag(hdr, line, field, nullptr, 0) == 1;
}

/**
 * @brief Read the GT field
 */
void getGT(bcf_hdr_t *header, bcf1_t *line, int isample, int *gt, int &ngt, bool &phased)
{
    bcf_unpack(line, BCF_UN_FMT);
    bcf_fmt_t *fmt_ptr = bcf_get_fmt(header, line, "GT");
    if(fmt_ptr)
    {
        ngt = fmt_ptr->n;
        if(ngt > MAX_GT)
        {
            std::cerr << "[W] Found a variant with more " << ngt << " > " << MAX_GT
                      << " (max) alt alleles. These become no-calls."
                      << "\n";
            for(int i = 0; i < MAX_GT; ++i)
            {
                gt[i] = -1;
            }
            phased = false;
            return;
        }

        // check if all alleles are populated
        switch(fmt_ptr->type)
        {
            case BCF_BT_INT8:
            {
                static const bcfhelpers::_impl::bcf_get_gts<int8_t> b;
                b(fmt_ptr, isample, gt, ngt, phased);
            }
                break;
            case BCF_BT_INT16:
            {
                static const bcfhelpers::_impl::bcf_get_gts<int16_t> b;
                b(fmt_ptr, isample, gt, ngt, phased);
            }
                break;
            case BCF_BT_INT32:
            {
                static const bcfhelpers::_impl::bcf_get_gts<int32_t> b;
                b(fmt_ptr, isample, gt, ngt, phased);
            }
                break;
            default:
                error("Unsupported GT type: %d at %s:%d\n", fmt_ptr->type, header->id[BCF_DT_CTG][line->rid].key,
                      line->pos + 1);
                break;
        }
    }
    else
    {
        ngt = 0;
        phased = false;
    }
}

/** read GQ(X) -- will use in this order: GQX, GQ, -1
 *
 * This is somewhat Illumina-specific, TODO: refactor to avoid using this
 * function in a general setting
 * */
void getGQ(const bcf_hdr_t *header, bcf1_t *line, int isample, float &gq)
{
    using namespace _impl;
    static const bcf_get_numeric_format<float> gf;

    std::vector<float> values;
    gf(header, line, "GQX", isample, values);
    if(values.empty())
    {
        gf(header, line, "GQ", isample, values);
    }
    if(values.size() > 1)
    {
        std::cerr << "[W] too many GQ fields at "
                  << header->id[BCF_DT_CTG][line->rid].key << ":" << line->pos
                  << "\n";
    }
    if(values.empty())
    {
        gq = -1;
    }
    else
    {
        gq = values[0];
    }
}

/** read AD */
void getAD(const bcf_hdr_t *header, bcf1_t *line, int isample, int *ad, int max_ad)
{
    using namespace _impl;
    static const bcf_get_numeric_format<int> gf;

    std::vector<int> values;
    gf(header, line, "AD", isample, values);
    if(max_ad < (int) values.size())
    {
        std::cerr << "[W] too many AD fields at "
                  << header->id[BCF_DT_CTG][line->rid].key << ":" << line->pos
                  << " max_ad = " << max_ad << " retrieved: " << values.size()
                  << "\n";
    }
    for(size_t q = 0; q < values.size(); ++q)
    {
        ad[q] = values[q];
    }
}

/** read DP(I) -- will use in this order: DP, DPI, -1 */
void getDP(const bcf_hdr_t *header, bcf1_t *line, int isample, int &dp)
{
    using namespace _impl;
    static const bcf_get_numeric_format<int> gf;

    std::vector<int> values;
    gf(header, line, "DP", isample, values);
    if(values.empty())
    {
        gf(header, line, "DPI", isample, values);
    }
    if(values.size() > 1)
    {
        std::cerr << "[W] too many DP fields at "
                  << header->id[BCF_DT_CTG][line->rid].key << ":" << line->pos
                  << "\n";
    }
    if(values.empty())
    {
        dp = 0;
    }
    else
    {
        dp = values[0];
    }
}

/** read a format field as a single int. */
int getFormatInt(const bcf_hdr_t *header, bcf1_t *line, const char *field, int isample, int defaultresult)
{
    using namespace _impl;
    static const bcf_get_numeric_format<int> gf;

    std::vector<int> values;
    gf(header, line, field, isample, values);
    if(values.size() > 1)
    {
        std::ostringstream os;
        os << "[W] too many " << field << " fields at "
           << header->id[BCF_DT_CTG][line->rid].key << ":" << line->pos;
        throw importexception(os.str());
    }
    if(values.size() == 1)
    {
        return values[0];
    }
    return defaultresult;
}

std::vector<int> getFormatInts(const bcf_hdr_t *header, bcf1_t *line, const char *field, int isample)
{
    using namespace _impl;
    static const bcf_get_numeric_format<int> gf;
    std::vector<int> result;
    gf(header, line, field, isample, result);
    return result;
}

/** read a format field as a single float. default return value is NaN */
float getFormatFloat(const bcf_hdr_t *header, bcf1_t *line, const char *field, int isample)
{
    using namespace _impl;
    float result = std::numeric_limits<float>::quiet_NaN();
    static const bcf_get_numeric_format<float> gf;

    std::vector<float> values;
    gf(header, line, field, isample, values);
    if(values.size() > 1)
    {
        std::ostringstream os;
        os << "[W] too many " << field << " fields at "
           << header->id[BCF_DT_CTG][line->rid].key << ":" << line->pos;
        throw importexception(os.str());
    }
    if(values.size() == 1)
    {
        result = values[0];
    }
    return result;
}

std::vector<float> getFormatFloats(const bcf_hdr_t *header, bcf1_t *line, const char *field, int isample)
{
    using namespace _impl;
    static const bcf_get_numeric_format<float> gf;
    std::vector<float> result;
    gf(header, line, field, isample, result);
    return result;
}

/** read a format field as a single double. result will not be overwritten on
 * failure */
std::string getFormatString(const bcf_hdr_t *hdr, bcf1_t *line, const char *field, int isample, const char *result)
{
    int nsmpl = bcf_hdr_nsamples(hdr);
    if(isample >= nsmpl)
    {
        return result;
    }

    bcf_unpack(line, BCF_UN_FMT);
    int tag_id = bcf_hdr_id2int(hdr, BCF_DT_ID, field);

    if(!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, tag_id))
    {
        return result;
    }

    if(!(line->unpacked & BCF_UN_FMT))
    {
        bcf_unpack(line, BCF_UN_FMT);
    }

    // index in format fields
    int i = 0;
    for(i = 0; i < line->n_fmt; i++)
    {
        if(line->d.fmt[i].id == tag_id)
        {
            break;
        }
    }

    if(i == line->n_fmt)
    {
        return result;
    }

    bcf_fmt_t *fmt = &line->d.fmt[i];

    if(fmt == NULL)
    {
        return result;
    }

    int type = fmt->type;

    if(fmt->n < 1)
    {
        return result;
    }

    std::string str_result = result;
    if(type == BCF_BT_FLOAT)
    {
        static const float bcf_missing_float = missing_float();
        // cppcheck-suppress invalidPointerCast
        float res = *((float *) (fmt->p + isample * fmt->size));
        if(res != bcf_missing_float)
        {
            str_result = std::to_string(res);
        }
    }
    else if(type == BCF_BT_INT8)
    {
        int8_t r = *((int8_t *) (fmt->p + isample * fmt->size));
        if(r != bcf_int8_missing && r != bcf_int8_vector_end)
        {
            str_result = std::to_string(r);
        }
    }
    else if(type == BCF_BT_INT16)
    {
        int16_t r = *((int16_t *) (fmt->p + isample * fmt->size));
        if(r != bcf_int16_missing && r != bcf_int16_vector_end)
        {
            str_result = std::to_string(r);
        }
    }
    else if(type == BCF_BT_INT32)
    {
        int32_t r = *((int32_t *) (fmt->p + isample * fmt->size));
        if(r != bcf_int32_missing && r != bcf_int32_vector_end)
        {
            str_result = std::to_string(r);
        }
    }
    else
    {
        const char *src = (const char *) fmt->p + isample * fmt->size;
        if(src && fmt->size)
        {
            str_result = std::string(src, (unsigned long) fmt->size);
            // deal with 0 padding
            str_result.resize(strlen(str_result.c_str()));
        }
    }

    return str_result;
}

/** update format string for a single sample.  */
void setFormatStrings(const bcf_hdr_t *hdr, bcf1_t *line, const char *field, const std::vector<std::string> &formats)
{
    // TODO this can probably be done faster / better
    std::unique_ptr<const char *[]> p_fmts = std::unique_ptr<const char *[]>(new const char *[line->n_sample]);
    bcf_unpack(line, BCF_UN_FMT);
    bool any_nonempty = false;
    if(formats.size() != line->n_sample)
    {
        std::ostringstream os;
        os << "[W] cannot update format " << field << " "
           << hdr->id[BCF_DT_ID][line->rid].key << ":" << line->pos
           << " -- we have " << formats.size() << " for " << line->n_sample
           << " samples";
        throw importexception(os.str());
    }
    for(int si = 0; si < line->n_sample; ++si)
    {
        p_fmts.get()[si] = formats[si].c_str();
        if(!formats[si].empty())
        {
            any_nonempty = true;
        }
    }
    int res;
    if(any_nonempty)
    {
        res = bcf_update_format_string(hdr, line, field, p_fmts.get(), line->n_sample);
    }
    else
    {
        res = bcf_update_format_string(hdr, line, field, NULL, 0);
    }

    if(res != 0)
    {
        std::ostringstream os;
        os << "[W] cannot update format " << field << " "
           << hdr->id[BCF_DT_ID][line->rid].key << ":" << line->pos
           << " -- we have " << formats.size() << " for " << line->n_sample
           << " samples";
        throw importexception(os.str());
    }
}

/** update format with single float values.  */
void setFormatFloats(const bcf_hdr_t *header, bcf1_t *line, const char *field, const std::vector<float> &value)
{
    if(value.empty())
    {
        int res = bcf_update_format(header, line, field, NULL, 0, 0);
        if(res != 0)
        {
            std::ostringstream os;
            os << "[W] cannot update format " << field << " "
               << header->id[BCF_DT_CTG][line->rid].key << ":" << line->pos
               << " -- we have " << value.size() << " for " << line->n_sample
               << " samples";

            throw importexception(os.str());
        }
        return;
    }
    std::unique_ptr<float[]> p_dbl = std::unique_ptr<float[]>(new float[line->n_sample]);

    for(size_t i = 0; i < line->n_sample; ++i)
    {
        if(i < value.size())
        {
            p_dbl.get()[i] = value[i];
        }
        else
        {
            p_dbl.get()[i] = std::numeric_limits<float>::quiet_NaN();
        }
    }

    int res = bcf_update_format_float(header, line, field, p_dbl.get(), line->n_sample);
    if(res != 0)
    {
        std::ostringstream os;
        os << "[W] cannot update format " << field << " "
           << header->id[BCF_DT_CTG][line->rid].key << ":" << line->pos
           << " -- we have " << value.size() << " for " << line->n_sample
           << " samples";
        throw importexception(os.str());
    }
}

void setFormatInts(const bcf_hdr_t *header, bcf1_t *line, const char *field, const std::vector<int> &value)
{
    if(value.empty())
    {
        int res = bcf_update_format(header, line, field, NULL, 0, 0);
        if(res != 0)
        {
            std::ostringstream os;
            os << "[W] cannot update format " << field << " "
               << header->id[BCF_DT_CTG][line->rid].key << ":" << line->pos
               << " -- we have " << value.size() << " for " << line->n_sample
               << " samples";

            throw importexception(os.str());
        }
        return;
    }
    std::unique_ptr<int[]> p_dbl = std::unique_ptr<int[]>(new int[line->n_sample]);

    for(size_t i = 0; i < line->n_sample; ++i)
    {
        if(i < value.size())
        {
            p_dbl.get()[i] = value[i];
        }
        else
        {
            p_dbl.get()[i] = bcf_int32_missing;
        }
    }

    int res = bcf_update_format_int32(header, line, field, p_dbl.get(), line->n_sample);
    if(res != 0)
    {
        std::ostringstream os;
        os << "[W] cannot update format " << field << " "
           << header->id[BCF_DT_CTG][line->rid].key << ":" << line->pos
           << " -- we have " << value.size() << " for " << line->n_sample
           << " samples";
        throw importexception(os.str());
    }
}

/** return sample names from header */
std::list<std::string> getSampleNames(const bcf_hdr_t *hdr)
{
    std::list<std::string> l;
    for(int i = 0; i < bcf_hdr_nsamples(hdr); ++i)
    {
        std::string samplename = hdr->samples[i];
        if(samplename == "*")
        {
            std::cerr << "Skipping sample named '*'"
                      << "\n";
            continue;
        }
        l.push_back(samplename);
    }
    return l;
}

/** return number of reference padding bases */
int addRefPad(bcf_hdr_t *hdr, bcf1_t *rec, FastaFile const &ref, int npad)
{
    if(npad <= 0)
    {
        error("npad <= 0");
    }
    int64_t start, end;
    getLocation(hdr, rec, start, end);
    std::string pad = ref.query(bcf_hdr_int2id(hdr, BCF_DT_CTG, rec->rid), start - npad, start - 1);
    rec->pos -= npad;
    std::string new_alleles = "";
    for(int i = 0; i < rec->n_allele; i++)
    {
        auto type = classifyAlleleString(rec->d.allele[i]).first;
        if(type == AlleleType::NUC)
        {
            if(i == 0)
            {
                new_alleles += pad;
                new_alleles += (std::string) rec->d.allele[i];
            }
            else
            {
                new_alleles += ",";
                new_alleles += pad;
                new_alleles += (std::string) rec->d.allele[i];
            }
        }
        else if(type == AlleleType::SYMBOLIC_DEL)
        {
            new_alleles += ",";
            new_alleles += pad;
        }
        else
        {
            new_alleles += ",";
            new_alleles += (std::string) rec->d.allele[i];
        }
    }
    bcf_update_alleles_str(hdr, rec, new_alleles.c_str());

    std::string cigar = getInfoString(hdr, rec, "CIGAR", "");
    if(!cigar.empty())
    {
        size_t pos = 0;
        int matches = npad;
        while(pos < cigar.size() && cigar[pos] >= '0' && cigar[pos] <= '9')
        {
            ++pos;
        }
        if(pos > 0 && pos < cigar.size() && cigar[pos] == 'M')
        {
            matches += atoi(cigar.substr(0, pos).c_str());
            cigar = cigar.substr(pos + 1);
        }
        cigar = std::to_string(matches) + "M" + cigar;
        bcf_update_info_string(hdr, rec, "CIGAR", cigar.c_str());
    }

    return (npad);
}

/** return number of reference padding bases */
int isRefPadded(bcf1_t *line)
{
    bcf_unpack(line, BCF_UN_SHR);
    if(line->n_allele == 1)
    {
        return 0;
    }

    const char *ref = line->d.allele[0];
    const int reflen = (int) strlen(ref);

    int max_match = reflen;
    for(int al = 1; al < line->n_allele; ++al)
    {
        const char *alt = line->d.allele[al];
        auto type = classifyAlleleString(alt).first;
        // overlapping other variants don't affect ref padding
        if(type == AlleleType::OVERLAPPING_DEL || type == AlleleType::MARGINALIZED)
        {
            continue;
        }
        if(type != AlleleType::NUC)
        {
            max_match = 0;
            break;
        }
        int rpos = 0;
        for(rpos = 0; rpos < reflen; ++rpos, ++alt)
        {
            const char rb = *(ref + rpos);
            if(*alt == 0 || *alt != rb)
            {
                break;
            }
        }
        max_match = std::min(rpos, max_match);
    }
    return max_match;
}

// everything after here is modified from bcftools 1.3.1. vcfnorm.c

static void
split_info_numeric(const bcf_hdr_t *hdr, bcf1_t *src, bcf_info_t *info, int ialt, bcf1_t *dst, uint8_t *&tmp_arr1,
                   int &ntmp_arr1)
{
#define BRANCH_NUMERIC(type, type_t)                                    \
    {                                                                   \
        const char* tag = bcf_hdr_int2id(hdr, BCF_DT_ID, info->key);    \
        int ntmp = ntmp_arr1 / sizeof(type_t);                          \
        int ret = bcf_get_info_##type(hdr, src, tag, &tmp_arr1, &ntmp); \
        ntmp_arr1 = ntmp * sizeof(type_t);                              \
        assert(ret > 0);                                                \
        type_t* vals = (type_t*)tmp_arr1;                               \
        int len = bcf_hdr_id2length(hdr, BCF_HL_INFO, info->key);       \
        if (len == BCF_VL_A)                                            \
        {                                                               \
            assert(ret == src->n_allele - 1);                           \
            bcf_update_info_##type(hdr, dst, tag, vals + ialt, 1);      \
        }                                                               \
        else if (len == BCF_VL_R)                                       \
        {                                                               \
            assert(ret == src->n_allele);                               \
            if (ialt != 0)                                              \
                vals[1] = vals[ialt + 1];                               \
            bcf_update_info_##type(hdr, dst, tag, vals, 2);             \
        }                                                               \
        else if (len == BCF_VL_G)                                       \
        {                                                               \
            assert(ret == src->n_allele * (src->n_allele + 1) / 2);     \
            if (ialt != 0)                                              \
            {                                                           \
                vals[1] = vals[bcf_alleles2gt(0, ialt + 1)];            \
                vals[2] = vals[bcf_alleles2gt(ialt + 1, ialt + 1)];     \
            }                                                           \
            bcf_update_info_##type(hdr, dst, tag, vals, 3);             \
        }                                                               \
        else                                                            \
            bcf_update_info_##type(hdr, dst, tag, vals, ret);           \
    }
    switch(bcf_hdr_id2type(hdr, BCF_HL_INFO, info->key))
    {
        case BCF_HT_INT:
        BRANCH_NUMERIC(int32, int32_t);
            break;
            // cppcheck-suppress invalidPointerCast
        case BCF_HT_REAL:
        BRANCH_NUMERIC(float, float);
            break;
    }
#undef BRANCH_NUMERIC
}
// Find n-th field in a comma-separated list and move it to dst.
// The memory areas may overlap.
#define STR_MOVE_NTH(dst, src, end, nth, len) \
    {                                         \
        char *ss = src, *se = src;            \
        int j = 0;                            \
        while (*se && se < (end))             \
        {                                     \
            if (*se == ',')                   \
            {                                 \
                if (j == nth)                 \
                    break;                    \
                j++;                          \
                ss = se + 1;                  \
            }                                 \
            se++;                             \
        }                                     \
        if (j == nth)                         \
        {                                     \
            int n = se - ss;                  \
            memmove((dst), ss, n);            \
            src = se;                         \
            len += n;                         \
        }                                     \
        else                                  \
            len = -1;                         \
    }

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wshadow"

static void
split_info_string(const bcf_hdr_t *hdr, bcf1_t *src, bcf_info_t *info, int ialt, bcf1_t *dst, uint8_t *&tmp_arr1,
                  int &ntmp_arr1)
{
    const char *tag = bcf_hdr_int2id(hdr, BCF_DT_ID, info->key);
    int ret = bcf_get_info_string(hdr, src, tag, &tmp_arr1, &ntmp_arr1);
    assert(ret > 0);

    kstring_t str;
    str.m = (size_t) ntmp_arr1;
    str.l = (size_t) ret;
    str.s = (char *) tmp_arr1;

    int len = bcf_hdr_id2length(hdr, BCF_HL_INFO, info->key);
    if(len == BCF_VL_A)
    {
        char *tmp = str.s;
        int len = 0;
        STR_MOVE_NTH(str.s, tmp, str.s + str.l, ialt, len);
        if(len < 0)
        {
            return;
        } // wrong number of fields: skip
        str.s[len] = 0;
        bcf_update_info_string(hdr, dst, tag, str.s);
    }
    else if(len == BCF_VL_R)
    {
        char *tmp = str.s;
        int len = 0;
        STR_MOVE_NTH(str.s, tmp, str.s + str.l, 0, len);
        str.s[len] = ',';
        tmp++;
        len++;
        STR_MOVE_NTH(&str.s[len], tmp, str.s + str.l, ialt, len);
        if(len < 0)
        {
            return;
        } // wrong number of fields: skip
        str.s[len] = 0;
        bcf_update_info_string(hdr, dst, tag, str.s);
    }
    else if(len == BCF_VL_G)
    {
        // cppcheck-suppress duplicateExpression
        int i0a = bcf_alleles2gt(0, ialt + 1),
                iaa = bcf_alleles2gt(ialt + 1, ialt + 1);
        char *tmp = str.s;
        int len = 0;
        STR_MOVE_NTH(str.s, tmp, str.s + str.l, 0, len);
        str.s[len] = ',';
        tmp++;
        len++;
        STR_MOVE_NTH(&str.s[len], tmp, str.s + str.l, i0a - 1, len);
        if(len < 0)
        {
            return;
        } // wrong number of fields: skip
        str.s[len] = ',';
        tmp++;
        len++;
        STR_MOVE_NTH(&str.s[len], tmp, str.s + str.l, iaa - i0a - 1, len);
        if(len < 0)
        {
            return;
        } // wrong number of fields: skip
        str.s[len] = 0;
        bcf_update_info_string(hdr, dst, tag, str.s);
    }
    else
        bcf_update_info_string(hdr, dst, tag, str.s);
}

#pragma clang diagnostic pop

static void
split_info_flag(const bcf_hdr_t *hdr, bcf1_t *src, bcf_info_t *info, int ialt, bcf1_t *dst, uint8_t *&tmp_arr1,
                int &ntmp_arr1)
{
    const char *tag = bcf_hdr_int2id(hdr, BCF_DT_ID, info->key);
    int ret = bcf_get_info_flag(hdr, src, tag, &tmp_arr1, &ntmp_arr1);
    bcf_update_info_flag(hdr, dst, tag, NULL, ret);
}

static void
split_format_genotype(const bcf_hdr_t *hdr, bcf1_t *src, bcf_fmt_t *fmt, int ialt, bcf1_t *dst, uint8_t *&tmp_arr1,
                      int &ntmp_arr1)
{
    int ntmp = ntmp_arr1 / 4;
    int ngts = bcf_get_genotypes(hdr, src, &tmp_arr1, &ntmp);
    ntmp_arr1 = ntmp * 4;
    assert(ngts > 0);

    int32_t *gt = (int32_t *) tmp_arr1;
    int i, j, nsmpl = bcf_hdr_nsamples(hdr);
    ngts /= nsmpl;
    for(i = 0; i < nsmpl; i++)
    {
        for(j = 0; j < ngts; j++)
        {
            if(gt[j] == bcf_int32_vector_end)
            {
                break;
            }
            if(bcf_gt_is_missing(gt[j]) || bcf_gt_allele(gt[j]) == 0)
            {
                continue; // missing allele or ref: leave as is
            }
            if(bcf_gt_allele(gt[j]) == ialt + 1)
            {
                gt[j] = bcf_gt_unphased(1) | bcf_gt_is_phased(gt[j]); // set to first ALT
            }
            else
            {
                gt[j] = bcf_gt_unphased(0) | bcf_gt_is_phased(gt[j]); // set to REF
            }
        }
        gt += ngts;
    }
    bcf_update_genotypes(hdr, dst, tmp_arr1, ngts * nsmpl);
}

static void
split_format_numeric(const bcf_hdr_t *hdr, bcf1_t *src, bcf_fmt_t *fmt, int ialt, bcf1_t *dst, uint8_t *&tmp_arr1,
                     int &ntmp_arr1)
{
#define BRANCH_NUMERIC(type, type_t, is_vector_end, set_vector_end)                                                                                             \
    {                                                                                                                                                           \
        const char* tag = bcf_hdr_int2id(hdr, BCF_DT_ID, fmt->id);                                                                                              \
        int ntmp = ntmp_arr1 / sizeof(type_t);                                                                                                                  \
        int nvals = bcf_get_format_##type(hdr, src, tag, &tmp_arr1, &ntmp);                                                                                     \
        ntmp_arr1 = ntmp * sizeof(type_t);                                                                                                                      \
        assert(nvals > 0);                                                                                                                                      \
        type_t* vals = (type_t*)tmp_arr1;                                                                                                                       \
        int len = bcf_hdr_id2length(hdr, BCF_HL_FMT, fmt->id);                                                                                                  \
        int i, nsmpl = bcf_hdr_nsamples(hdr);                                                                                                                   \
        if (nvals == nsmpl) /* all values are missing */                                                                                                        \
        {                                                                                                                                                       \
            bcf_update_format_##type(hdr, dst, tag, vals, nsmpl);                                                                                               \
            return;                                                                                                                                             \
        }                                                                                                                                                       \
        if (len == BCF_VL_A)                                                                                                                                    \
        {                                                                                                                                                       \
            assert(nvals == (src->n_allele - 1) * nsmpl);                                                                                                       \
            nvals /= nsmpl;                                                                                                                                     \
            type_t *src_vals = vals, *dst_vals = vals;                                                                                                          \
            for (i = 0; i < nsmpl; i++)                                                                                                                         \
            {                                                                                                                                                   \
                dst_vals[0] = src_vals[ialt];                                                                                                                   \
                dst_vals += 1;                                                                                                                                  \
                src_vals += nvals;                                                                                                                              \
            }                                                                                                                                                   \
            bcf_update_format_##type(hdr, dst, tag, vals, nsmpl);                                                                                               \
        }                                                                                                                                                       \
        else if (len == BCF_VL_R)                                                                                                                               \
        {                                                                                                                                                       \
            assert(nvals == src->n_allele * nsmpl);                                                                                                             \
            nvals /= nsmpl;                                                                                                                                     \
            type_t *src_vals = vals, *dst_vals = vals;                                                                                                          \
            for (i = 0; i < nsmpl; i++)                                                                                                                         \
            {                                                                                                                                                   \
                dst_vals[0] = src_vals[0];                                                                                                                      \
                dst_vals[1] = src_vals[ialt + 1];                                                                                                               \
                dst_vals += 2;                                                                                                                                  \
                src_vals += nvals;                                                                                                                              \
            }                                                                                                                                                   \
            bcf_update_format_##type(hdr, dst, tag, vals, nsmpl * 2);                                                                                           \
        }                                                                                                                                                       \
        else if (len == BCF_VL_G)                                                                                                                               \
        {                                                                                                                                                       \
            if (nvals != src->n_allele * (src->n_allele + 1) / 2 * nsmpl && nvals != src->n_allele * nsmpl)                                                     \
                error("Error at %s:%d, the tag %s has wrong number of fields\n", bcf_seqname(hdr, src), src->pos + 1, bcf_hdr_int2id(hdr, BCF_DT_ID, fmt->id)); \
            nvals /= nsmpl;                                                                                                                                     \
            int all_haploid = nvals == src->n_allele ? 1 : 0;                                                                                                   \
            type_t *src_vals = vals, *dst_vals = vals;                                                                                                          \
            for (i = 0; i < nsmpl; i++)                                                                                                                         \
            {                                                                                                                                                   \
                int haploid = all_haploid;                                                                                                                      \
                if (!haploid)                                                                                                                                   \
                {                                                                                                                                               \
                    int j;                                                                                                                                      \
                    for (j = 0; j < nvals; j++)                                                                                                                 \
                        if (is_vector_end)                                                                                                                      \
                            break;                                                                                                                              \
                    if (j != nvals)                                                                                                                             \
                        haploid = 1;                                                                                                                            \
                }                                                                                                                                               \
                dst_vals[0] = src_vals[0];                                                                                                                      \
                if (haploid)                                                                                                                                    \
                {                                                                                                                                               \
                    dst_vals[1] = src_vals[ialt + 1];                                                                                                           \
                    if (!all_haploid)                                                                                                                           \
                        set_vector_end;                                                                                                                         \
                }                                                                                                                                               \
                else                                                                                                                                            \
                {                                                                                                                                               \
                    dst_vals[1] = src_vals[bcf_alleles2gt(0, ialt + 1)];                                                                                        \
                    dst_vals[2] = src_vals[bcf_alleles2gt(ialt + 1, ialt + 1)];                                                                                 \
                }                                                                                                                                               \
                dst_vals += all_haploid ? 2 : 3;                                                                                                                \
                src_vals += nvals;                                                                                                                              \
            }                                                                                                                                                   \
            bcf_update_format_##type(hdr, dst, tag, vals, all_haploid ? nsmpl * 2 : nsmpl * 3);                                                                 \
        }                                                                                                                                                       \
        else                                                                                                                                                    \
            bcf_update_format_##type(hdr, dst, tag, vals, nvals);                                                                                               \
    }
    switch(bcf_hdr_id2type(hdr, BCF_HL_FMT, fmt->id))
    {
        case BCF_HT_INT:
        BRANCH_NUMERIC(int32, int32_t, src_vals[j] == bcf_int32_vector_end, dst_vals[2] = bcf_int32_vector_end);
            break;
        case BCF_HT_REAL:
        BRANCH_NUMERIC(float, float, bcf_float_is_vector_end(src_vals[j]), bcf_float_set_vector_end(dst_vals[2]));
            break;
    }
#undef BRANCH_NUMERIC
}

static void squeeze_format_char(char *str, int src_blen, int dst_blen, int n)
{
    int i, isrc = 0, idst = 0;
    for(i = 0; i < n; i++)
    {
        memmove(str + idst, str + isrc, dst_blen);
        idst += dst_blen;
        isrc += src_blen;
    }
}

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wshadow"

static void
split_format_string(const bcf_hdr_t *hdr, bcf1_t *src, bcf_fmt_t *fmt, int ialt, bcf1_t *dst, uint8_t *&tmp_arr1,
                    int &ntmp_arr1)
{
    const char *tag = bcf_hdr_int2id(hdr, BCF_DT_ID, fmt->id);
    int ret = bcf_get_format_char(hdr, src, tag, &tmp_arr1, &ntmp_arr1);
    assert(ret > 0);

    kstring_t str;
    str.m = (size_t) ntmp_arr1;
    str.l = (size_t) ret;
    str.s = (char *) tmp_arr1;

    int nsmpl = bcf_hdr_nsamples(hdr);
    int len = bcf_hdr_id2length(hdr, BCF_HL_FMT, fmt->id);
    if(len == BCF_VL_A)
    {
        int i, blen = ret / nsmpl, maxlen = 0;
        char *ptr = str.s;
        for(i = 0; i < nsmpl; i++)
        {
            char *tmp = ptr;
            int len = 0;
            STR_MOVE_NTH(tmp, tmp, ptr + blen, ialt, len);
            if(len < 0)
            {
                return;
            } // wrong number of fields: skip
            if(maxlen < len)
            {
                maxlen = len;
            }
            ptr += blen;
        }
        if(maxlen < blen)
        {
            squeeze_format_char(str.s, blen, maxlen, nsmpl);
        }
        bcf_update_format_char(hdr, dst, tag, str.s, nsmpl * maxlen);
    }
    else if(len == BCF_VL_R)
    {
        int i, blen = ret / nsmpl, maxlen = 0;
        char *ptr = str.s;
        for(i = 0; i < nsmpl; i++)
        {
            char *tmp = ptr;
            int len = 0;
            STR_MOVE_NTH(ptr, tmp, ptr + blen, 0, len);
            ptr[len] = ',';
            tmp++;
            len++;
            STR_MOVE_NTH(&ptr[len], tmp, ptr + blen, ialt, len);
            if(len < 0)
            {
                return;
            } // wrong number of fields: skip
            if(maxlen < len)
            {
                maxlen = len;
            }
            ptr += blen;
        }
        if(maxlen < blen)
        {
            squeeze_format_char(str.s, blen, maxlen, nsmpl);
        }
        bcf_update_format_char(hdr, dst, tag, str.s, nsmpl * maxlen);
    }
    else if(len == BCF_VL_G)
    {
        // cppcheck-suppress duplicateExpression
        int i, blen = ret / nsmpl, maxlen = 0, i0a = bcf_alleles2gt(0, ialt + 1),
                iaa = bcf_alleles2gt(ialt + 1, ialt + 1);
        char *ptr = str.s;
        for(i = 0; i < nsmpl; i++)
        {
            char *se = ptr, *sx = ptr + blen;
            int nfields = 1;
            while(*se && se < sx)
            {
                if(*se == ',')
                {
                    nfields++;
                }
                se++;
            }
            assert(nfields == src->n_allele * (src->n_allele + 1) / 2 || nfields == src->n_allele);
            int len = 0;
            if(nfields == src->n_allele) // haploid
            {
                char *tmp = ptr;
                STR_MOVE_NTH(&ptr[len], tmp, ptr + blen, 0, len);
                ptr[len] = ',';
                tmp++;
                len++;
                STR_MOVE_NTH(&ptr[len], tmp, ptr + blen, ialt, len);
                if(len < 0)
                {
                    return;
                } // wrong number of fields: skip
            }
            else // diploid
            {
                char *tmp = ptr;
                STR_MOVE_NTH(&ptr[len], tmp, ptr + blen, 0, len);
                ptr[len] = ',';
                tmp++;
                len++;
                STR_MOVE_NTH(&ptr[len], tmp, ptr + blen, i0a - 1, len);
                if(len < 0)
                {
                    return;
                } // wrong number of fields: skip
                ptr[len] = ',';
                tmp++;
                len++;
                STR_MOVE_NTH(&ptr[len], tmp, ptr + blen, iaa - i0a - 1, len);
                if(len < 0)
                {
                    return;
                } // wrong number of fields: skip
            }
            if(maxlen < len)
            {
                maxlen = len;
            }
            ptr += blen;
        }
        if(maxlen < blen)
        {
            squeeze_format_char(str.s, blen, maxlen, nsmpl);
        }
        bcf_update_format_char(hdr, dst, tag, str.s, nsmpl * maxlen);
    }
    else
    {
        bcf_update_format_char(hdr, dst, tag, str.s, (int) str.l);
    }
}

#pragma clang diagnostic pop

void splitMultiAllelics(bcf_hdr_t *hdr, bcf1_t *line, std::vector<p_bcf1> &out)
{
    bcf_unpack(line, BCF_UN_ALL);

    // Init the target biallelic lines
    int ntmp_lines = line->n_allele - 1;
    kstring_t tmp = {0, 0, 0};
    kputs(line->d.allele[0], &tmp);
    kputc(',', &tmp);
    int rlen = (int) tmp.l;
    int gt_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "GT");

    // work arrays - this should probably be moved outside the function.
    uint8_t *tmp_arr1 = NULL;
    int ntmp_arr1 = 0;

    out.clear();
    for(int i = 0; i < ntmp_lines; i++) // for each ALT allele
    {
        p_bcf1 p_dst = p(bcf_init1());
        bcf1_t *dst = p_dst.get();
        out.push_back(p_dst);

        dst->rid = line->rid;
        dst->pos = line->pos;
        dst->qual = line->qual;

        // Not quite sure how to handle IDs, they can be assigned to a specific
        // ALT.  For now we leave the ID unchanged for all.
        bcf_update_id(hdr, dst, line->d.id ? line->d.id : ".");

        tmp.l = (size_t) rlen;
        kputs(line->d.allele[i + 1], &tmp);
        bcf_update_alleles_str(hdr, dst, tmp.s);

        if(line->d.n_flt)
        {
            bcf_update_filter(hdr, dst, line->d.flt, line->d.n_flt);
        }

        for(int j = 0; j < line->n_info; j++)
        {
            bcf_info_t *info = &line->d.info[j];
            int type = bcf_hdr_id2type(hdr, BCF_HL_INFO, info->key);
            if(type == BCF_HT_INT || type == BCF_HT_REAL)
            {
                split_info_numeric(hdr, line, info, i, dst, tmp_arr1, ntmp_arr1);
            }
            else if(type == BCF_HT_FLAG)
            {
                split_info_flag(hdr, line, info, i, dst, tmp_arr1, ntmp_arr1);
            }
            else
            {
                split_info_string(hdr, line, info, i, dst, tmp_arr1, ntmp_arr1);
            }
        }

        dst->n_sample = line->n_sample;
        for(int j = 0; j < line->n_fmt; j++)
        {
            bcf_fmt_t *fmt = &line->d.fmt[j];
            int type = bcf_hdr_id2type(hdr, BCF_HL_FMT, fmt->id);
            if(fmt->id == gt_id)
            {
                split_format_genotype(hdr, line, fmt, i, dst, tmp_arr1, ntmp_arr1);
            }
            else if(type == BCF_HT_INT || type == BCF_HT_REAL)
            {
                split_format_numeric(hdr, line, fmt, i, dst, tmp_arr1, ntmp_arr1);
            }
            else
            {
                split_format_string(hdr, line, fmt, i, dst, tmp_arr1, ntmp_arr1);
            }
        }
    }
    free(tmp.s);
    free(tmp_arr1);
}

void addSymbolic(const bcf_hdr_t *hdr, bcf1_t *rec)
{
    assert(rec->n_allele == 2);
    std::string new_alleles = (std::string) rec->d.allele[0] + "," + (std::string) rec->d.allele[1] + ",<M>";
    bcf_update_alleles_str(hdr, rec, new_alleles.c_str());
    bcf_unpack(rec, BCF_UN_ALL);
    for(int i = 0; i < rec->n_info; i++)
    {
        int nval = 2;
        bcf_info_t *info = &rec->d.info[i];
        int type = bcf_hdr_id2type(hdr, BCF_HL_INFO, info->key);
        const char *tag = bcf_hdr_int2id(hdr, BCF_DT_ID, info->key);
        int len = bcf_hdr_id2length(hdr, BCF_HL_INFO, info->key);
        if(len == BCF_VL_A)
        {
            if(type == BCF_HT_INT)
            {
                int32_t *new_val = (int32_t *) malloc(2 * sizeof(int32_t));
                bcf_get_info_int32(hdr, rec, tag, &new_val, &nval);
                new_val[1] = bcf_int32_missing;
                bcf_update_info_int32(hdr, rec, tag, new_val, 2);
                free(new_val);
            }
            else if(type == BCF_HT_REAL)
            {
                float *new_val = (float *) malloc(2 * sizeof(float));
                bcf_get_info_int32(hdr, rec, tag, &new_val, &nval);
                bcf_float_set_missing(new_val[1]);
                bcf_update_info_float(hdr, rec, tag, new_val, 2);
                free(new_val);
            }
            else if(type == BCF_HT_FLAG)
            {
                continue;
            }
            else // string
            {
                nval = 0;
                char *oldstring = NULL;
                bcf_get_info_string(hdr, rec, tag, &oldstring, &nval);
                std::string newstring = (std::string) (oldstring) + ",.";
                bcf_update_info_string(hdr, rec, tag, newstring.c_str());
                free(oldstring);
            }
        }
    }
}

std::pair<AlleleType, std::string>
classifyAlleleString(std::string allele_string)
{
    stringutil::toUpper(allele_string);

    if(allele_string == "<DEL>")
    {
        return std::make_pair(AlleleType::SYMBOLIC_DEL, std::string());
    }
    else if(allele_string == "*")
    {
        return std::make_pair(AlleleType::OVERLAPPING_DEL, std::string("*"));
    }
    else if(allele_string == "<M>" || allele_string == "<*>")
    {
        return std::make_pair(AlleleType::MARGINALIZED, std::string("<M>"));
    }
    else if(allele_string == ".")
    {
        return std::make_pair(AlleleType::MISSING, std::string("."));
    }
    else if(allele_string.size() > 3 &&
            (allele_string.find('<') == 0 || allele_string.find('>') == allele_string.size() - 1))
    {
        return std::make_pair(AlleleType::SYMBOLIC_OTHER, allele_string);
    }
    else
    {
        // check alleles
        for(size_t qq = 0; qq < allele_string.size(); ++qq)
        {
            /* A	A	Adenine */
            /* C	C	Cytosine */
            /* G	G	Guanine */
            /* T	T	Thymine */
            /* U	U	Uracil */
            /* R	A or G	puRine */
            /* Y	C, T or U	pYrimidines */
            /* K	G, T or U	bases which are Ketones */
            /* M	A or C	bases with aMino groups */
            /* S	C or G	Strong interaction */
            /* W	A, T or U	Weak interaction */
            /* B	not A (i.e. C, G, T or U)	B comes after A */
            /* D	not C (i.e. A, G, T or U)	D comes after C */
            /* H	not G (i.e., A, C, T or U)	H comes after G */
            /* V	neither T nor U (i.e. A, C or G)	V comes after U */
            /* N	A C G T U	Nucleic acid */
            /* X	masked */
            /* -/.	gap of indeterminate length - not supported */
            if(allele_string[qq] != 'A' && allele_string[qq] != 'C' && allele_string[qq] != 'G' &&
               allele_string[qq] != 'T' && allele_string[qq] != 'U' && allele_string[qq] != 'R' &&
               allele_string[qq] != 'Y' && allele_string[qq] != 'K' && allele_string[qq] != 'M' &&
               allele_string[qq] != 'S' && allele_string[qq] != 'W' && allele_string[qq] != 'B' &&
               allele_string[qq] != 'D' && allele_string[qq] != 'H' && allele_string[qq] != 'V' &&
               allele_string[qq] != 'N' && allele_string[qq] != 'X')
            {
                return std::make_pair(AlleleType::UNKNOWN, allele_string);
            }
        }
        return std::make_pair(AlleleType::NUC, allele_string);
    }
}

std::ostream &operator<<(std::ostream &o, const AlleleType at)
{
    switch(at)
    {
        case AlleleType::NUC:
            o << "NUC";
            break;
        case AlleleType::MISSING:
            o << "MISSING";
            break;
        case AlleleType::SYMBOLIC_DEL:
            o << "SYMBOLIC_DEL";
            break;
        case AlleleType::OVERLAPPING_DEL:
            o << "OVERLAPPING_DEL";
            break;
        case AlleleType::MARGINALIZED:
            o << "MARGINALIZED";
            break;
        case AlleleType::SYMBOLIC_OTHER:
            o << "SYMBOLIC_OTHER";
            break;
        case AlleleType::UNKNOWN:
            o << "UNKNOWN";
            break;
    }
    return o;
}

/** check if a variant is a SNP */
bool isSNP(bcf_hdr_t *hdr, bcf1_t *rec)
{
    bcf_unpack(rec, BCF_UN_SHR);
    if(rec->n_allele < 2)
    {
        return false;
    }
    for(int i = 0; i < rec->n_allele; ++i)
    {
        auto type_and_str = classifyAlleleString(rec->d.allele[i]);
        if(i > 0 &&
           (type_and_str.first == AlleleType::OVERLAPPING_DEL || type_and_str.first == AlleleType::MARGINALIZED))
        {
            continue;
        }
        if(type_and_str.first != AlleleType::NUC || type_and_str.second.size() != 1)
        {
            return false;
        }
    }
    return true;
}

/**
 * Check if variant is a single-allelic MNP ie. can be trivially decomposed
 */
bool isMNP(bcf_hdr_t *hdr, bcf1_t *rec)
{
    bcf_unpack(rec, BCF_UN_SHR);
    if(rec->n_allele < 2)
    {
        return false;
    }
    size_t reflen = 0;
    for(int i = 0; i < rec->n_allele; i++)
    {
        auto type_and_str = classifyAlleleString(rec->d.allele[i]);
        if(i == 0 && type_and_str.first == AlleleType::NUC && type_and_str.second.size() > 1)
        {
            reflen = type_and_str.second.size();
        }
        if(i > 0 &&
           (type_and_str.first == AlleleType::OVERLAPPING_DEL || type_and_str.first == AlleleType::MARGINALIZED))
        {
            continue;
        }
        if(type_and_str.first != AlleleType::NUC || type_and_str.second.size() != reflen)
        {
            return false;
        }
    }
    return (true);
}

/**
 * Check if variant has any complex alleles which are neither insertions nor deletions
 *
 * This function detects reference padding
 */
bool isComplex(bcf_hdr_t *hdr, bcf1_t *rec)
{
    bcf_unpack(rec, BCF_UN_SHR);
    if(rec->n_allele <= 1)
    {
        return false;
    }
    if(isMNP(hdr, rec))
    {
        return true;
    }

    for(int i = 0; i < rec->n_allele; i++)
    {
        const char *ref = rec->d.allele[0];
        const char *alt = rec->d.allele[i];
        auto type = classifyAlleleString(alt).first;

        if(type != AlleleType::NUC
           && type != AlleleType::MARGINALIZED
           && type != AlleleType::SYMBOLIC_DEL
           && type != AlleleType::OVERLAPPING_DEL)
        {
            return true;
        }

        if(type != AlleleType::NUC || i == 0)
        {
            // skip overlapping DELs and REF allele
            continue;
        }

        while(*ref && *alt && *ref == *alt)
        {
            ++ref;
            ++alt;
        }

        // is it a deletion or insertion?
        if(!((strlen(ref) > 0 && strlen(alt) == 0) || (strlen(ref) == 0 && strlen(alt) > 0)))
        {
            return true;
        }
    }
    return false;
}

/**
 * Check if a variant has symbolic alleles beyond deletions and overlapping
 * or marginalized alleles
 */
bool isSymbolic(bcf_hdr_t *hdr, bcf1_t *rec)
{
    bcf_unpack(rec, BCF_UN_SHR);

    if(rec->n_allele <= 1)
    {
        return false;
    }

    for(int i = 0; i < rec->n_allele; i++)
    {
        const char *alt = rec->d.allele[i];
        auto type = classifyAlleleString(alt).first;

        if(type != AlleleType::NUC
           && type != AlleleType::MARGINALIZED
           && type != AlleleType::SYMBOLIC_DEL
           && type != AlleleType::OVERLAPPING_DEL)
        {
            return true;
        }
    }
    return false;
}

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
                 bool marginalize)
{
    assert((keep_alleles & ref_alleles) == 0);
    std::vector<int> old_ad = getFormatInts(hdr, in, "AD", isample);
    assert(in->n_allele <= 1 + 64 - std::min(clz64(keep_alleles), clz64(ref_alleles)));
    assert(old_ad.size() <= (size_t) (1 + 64 - std::min(clz64(keep_alleles), clz64(ref_alleles))));
    static const int zero = 0;
    new_ad.resize(popcount64(keep_alleles) + 1 + (marginalize ? 1 : 0), zero);
    int new_al = 1;
    int old_al = 1;
    new_ad[0] = old_ad[0];
    while(keep_alleles || ref_alleles)
    {
        if(old_al >= (int) old_ad.size())
        {
            break;
        }
        if(keep_alleles & 1)
        {
            new_ad[new_al] += old_ad[old_al];
            ++new_al;
        }
        else if(ref_alleles & 1)
        {
            new_ad[0] += old_ad[old_al];
        }
        else if(marginalize)
        {
            new_ad.back() += old_ad[old_al];
        }

        keep_alleles >>= 1;
        ref_alleles >>= 1;
        ++old_al;
    }
}

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
 *
 * @return a vector of lists which give the GT alleles in sorted order
 */
std::vector<std::list<int> > const &makePLAlleleLists(int ngt, int nal)
{
    static std::unordered_map<uint64_t, std::vector<std::list<int> > > precomputed_masks;

    const uint64_t id = ((uint64_t) ngt) + (((uint64_t) nal) << 32);

    auto it = precomputed_masks.find(id);
    if(it == precomputed_masks.end())
    {
        // create allele masks (see VCF SPEC)
        std::vector<std::list<int> > allele_masks;
        const std::function<void(int, int, std::list<int>)> makeMasks =
                [&makeMasks, &allele_masks](int p, int n, std::list<int> suffix)
                {
                    for(int a = 0; a <= n; ++a)
                    {
                        if(p == 1)
                        {
                            auto new_suffix = suffix;
                            new_suffix.push_front(a);
                            allele_masks.push_back(new_suffix);
                        }
                        else if(p > 1)
                        {
                            auto new_suffix = suffix;
                            new_suffix.push_front(a);
                            makeMasks(p - 1, a, new_suffix);
                        }
                    }
                };
        makeMasks(ngt, nal - 1, std::list<int>{});
        it = precomputed_masks.emplace(std::make_pair(id, allele_masks)).first;
    }
    return it->second;
}

/**
 * Retrieve GL values from a record, from PL if available
 * @param hdr BCF header
 * @param in record to read from
 * @param gl target array
 */
void getGL(bcf_hdr_t *hdr,
           bcf1_t *in,
           std::vector<float> &gl)
{
    int ngl, ngl_arr = 0;
    float *gls = NULL;

    ngl = bcf_get_format_float(hdr, in, "GL", &gls, &ngl_arr);

    if(ngl <= 0)
    {
        int npl, npl_arr = 0, *pls = NULL;
        npl = bcf_get_format_int32(hdr, in, "PL", &pls, &npl_arr);
        if(npl <= 0)
        {
            free(gls);
            gl.clear();
            return;
        }

        gls = (float *) malloc(sizeof(float) * npl);
        ngl = npl;

        for(int i = 0; i < ngl; i++)
        {
            gls[i] = powf(10.f, -pls[i] / 10.f);
        }

        free(pls);
    }

    int *gt_arr = NULL, ngt, ngt_arr = 0;
    ngt = bcf_get_genotypes(hdr, in, &gt_arr, &ngt_arr);
    free(gt_arr);
    const int ploidy = ngt / bcf_hdr_nsamples(hdr);
    auto const &masks = makePLAlleleLists(ploidy, in->n_allele);

    // check the number is correct
    if(masks.size() * bcf_hdr_nsamples(hdr) != (size_t) ngl)
    {
        free(gls);
        gl.clear();
        return;
    }

    gl.resize((unsigned long) ngl);
    for(int i = 0; i < ngl; i++)
    {
        gl[i] = gls[i];
    }
    free(gls);
}

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
                 bool marginalize)
{
    assert((keep_alleles & ref_alleles) == 0);

#ifdef DEBUG_PL_DECOMPOSE
    std::cerr << "keep " << std::hex << keep_alleles << " ref " << ref_alleles << std::endl;
#endif

    const int nal = in->n_allele;
    assert(nal <= 1 + 64 - std::min(clz64(keep_alleles), clz64(ref_alleles)));

    std::vector<float> old_gl = getFormatFloats(hdr, in, "GL", isample);

    if(old_gl.size() == 0)
    {
        std::vector<int> old_pl = getFormatInts(hdr, in, "PL", isample);
        for(size_t i = 0; i < old_pl.size(); i++)
        {
            const float tmp = powf(10.f, -old_pl[i] / 10.f);
            old_gl.push_back(tmp);
        }
    }

    int gt[MAX_GT], ngt;
    bool phased;
    getGT(hdr, in, isample, gt, ngt, phased);

    if(old_gl.size() == 0 || ngt == 0)
    {
        // no GL is present
        new_pl.clear();
        return;
    }

    size_t marginal_allele = 0;
    if(marginalize)
    {
        // TODO check if this works
        for(int j = 1; j < in->n_allele; ++j)
        {
            if(strcmp(in->d.allele[j], "<M>") == 0)
            {
                marginal_allele = (size_t) j;
                keep_alleles |= 1ull << (j - 1);
                break;
            }
        }

        if(marginal_allele == 0)
        {
            // e.g. 2 alleles kept out of 3:
            // alleles = A   C,G,T
            // keep      A   C,G
            // -> keep_alleles = 0b0011
            // -> popcount(keep_alleles) = 2
            // -> add new marginal allele:
            // -> A    C,G,<M>
            // -> index of marginal allele is 2 + 1 [ref] = 3
            marginal_allele = popcount64(keep_alleles) + 1;
            keep_alleles |= 1ull << ((int) marginal_allele - 1);
        }
    }

    auto old_masks = makePLAlleleLists(ngt, in->n_allele);
    assert(old_masks.size() == old_gl.size());

    auto new_masks = makePLAlleleLists(ngt, popcount64(keep_alleles) + 1);
    static const int zero = 0;
    new_pl.resize(new_masks.size(), zero);

    std::vector<uint64_t> old_to_new;
    old_to_new.resize(in->n_allele);
    int old_al = 1;
    size_t new_al = 1;
    old_to_new[0] = 0;
    while(keep_alleles || ref_alleles)
    {
        if(keep_alleles & 1)
        {
            old_to_new[old_al] = new_al++;
        }
        else if(ref_alleles & 1)
        {
            old_to_new[old_al] = 0;
        }
        else if(marginalize)
        {
            old_to_new[old_al] = marginal_allele;
        }
        else
        {
            error("unassigned alleles when decomposing");
        }
        ++old_al;
        keep_alleles >>= 1;
        ref_alleles >>= 1;
    }

    auto compare_lists = [](std::list<int> const &l1, std::list<int> const &l2) -> bool
    {
        if(l1.size() != l2.size())
        {
            return false;
        }
        else
        {
            auto i1 = l1.begin(), i2 = l2.begin();
            while(i1 != l1.end() && i2 != l2.end())
            {
                if(*i1 != *i2)
                {
                    return false;
                }
                ++i1;
                ++i2;
            }
            return true;
        }
    };

    auto sorted_ins = [](int value, std::list<int> &target)
    {
        auto it = target.begin();
        while(it != target.end() && value > *it)
        {
            ++it;
        }
        target.insert(it, value);
    };

    static const float feps = std::numeric_limits<float>::min();
    std::vector<float> new_gl;
    new_gl.resize(new_masks.size(), feps);
    float glsum = 0.f;

    for(size_t p = 0; p < old_masks.size(); ++p)
    {
#ifdef DEBUG_PL_DECOMPOSE
        std::cerr << " old GT ";
        for(auto old_allele : old_masks[p])
        {
            std::cerr << old_allele;
        }
        std::cerr << std::endl;
#endif

        std::list<int> new_alleles;
        for(auto old_allele : old_masks[p])
        {
            sorted_ins((int) old_to_new[old_allele], new_alleles);
        }

        int new_pl_idx = 0;
        for(auto new_alleles_matched : new_masks)
        {
            if(compare_lists(new_alleles, new_alleles_matched))
            {
#ifdef DEBUG_PL_DECOMPOSE
                std::cerr << "   ... matched new GT ";
                for(auto new_allele : new_alleles)
                {
                    std::cerr << new_allele;
                }
                std::cerr << std::endl;
#endif
                new_gl[new_pl_idx] += old_gl[p];
                glsum += old_gl[p];
            }
            ++new_pl_idx;
        }
    }

    for(size_t i = 0; i < new_gl.size(); i++)
    {
        const int tmp = (int) (-10.0f * log10f(new_gl[i] / glsum));
        new_pl[i] = tmp;
    }
}

void updateCigar(bcf_hdr_t *hdr, bcf1_t *rec)
{
    std::string cigar = getInfoString(hdr, rec, "CIGAR", "");
    if(!cigar.empty())
    {
        std::unique_ptr<Alignment> aln(makeAlignment("klibg"));
        aln->setQuery(rec->d.allele[0]);
        cigar = "";
        for(int i = 1; i < rec->n_allele; i++)
        {
            if(!cigar.empty())
            {
                cigar += ",";
            }

            auto type = classifyAlleleString(rec->d.allele[i]).first;
            if(type == AlleleType::NUC)
            {
                aln->setRef(rec->d.allele[i]);
                int i0, i1, e0, e1;
                std::string al_cigar;
                aln->getCigar(i0, i1, e0, e1, al_cigar);
                cigar += al_cigar;
            }
            else
            {
                cigar += ".";
            }
        }
        bcf_update_info_string(hdr, rec, "CIGAR", cigar.c_str());
    }
}

} // namespace bcfhelpers
