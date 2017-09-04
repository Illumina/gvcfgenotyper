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
 *  \brief Error handling helper
 *
 * \file Error.hh
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#pragma once

#include <cstdarg>
#include <cstdio>
#include <stdexcept>

#ifndef ERRMSG_MAX
#define ERRMSG_MAX 2048
#endif

#ifndef CYGWIN

#ifdef __GNUC__

#include <cstdio>
#include <cstdlib>
#include <execinfo.h>
#include <signal.h>

#endif

inline void error_stop(const char* f, int l, const char* format, ...)
{

#ifdef _DEBUG
    fprintf(stderr, "\n[DEBUG] Error at %s:%i\n", f, l);
    { // extra Debug print when errors occur
        va_list ap;
        va_start(ap, format);
        vfprintf(stderr, format, ap);
        va_end(ap);
        fprintf(stderr, "\n");
    }

// in GNUC, do backtrace
#ifdef __GNUC__
    void* array[20];

    size_t size;
    size = backtrace(array, 20);
    backtrace_symbols_fd(array, size, 2);
#endif
    fprintf(stderr, "\n[DEBUG END]\n");
#endif

    char errmsg[ERRMSG_MAX];
    va_list ap;
    va_start(ap, format);
    vsnprintf(errmsg, ERRMSG_MAX, format, ap);
    va_end(ap);
    throw std::runtime_error(errmsg);
}

#else // CYGWIN / windows version

#include <cstdio>
#include <cstdlib>
#include <imagehlp.h>
#include <inttypes.h>
#include <stdexcept>
#include <windows.h>

static inline int addr2line(void const* const addr)
{
    char addr2line_cmd[512] = { 0 };
    HMODULE hModule = GetModuleHandleW(NULL);
    WCHAR program_name[MAX_PATH];
    GetModuleFileNameW(hModule, program_name, MAX_PATH);

    /* have addr2line map the address to the relent line in the code */
    sprintf(addr2line_cmd, "addr2line -f -p -e `cygpath -u '%.256ls'` %p", program_name, addr);

    /* This will print a nicely formatted string specifying the
     function and source line of the address */
    return system(addr2line_cmd);
}

inline void error_stop(const char* f, int l, const char* format, ...)
{
#ifdef _DEBUG
    fprintf(stderr, "\n[DEBUG] Error at %s:%i\n", f, l);
    { // extra Debug print when errors occur
        va_list ap;
        va_start(ap, format);
        vfprintf(stderr, format, ap);
        va_end(ap);
        fprintf(stderr, "\n");
    }

    const int kmaxcallers = 62;

    void* callers[kmaxcallers];
    int count = CaptureStackBackTrace(0, kmaxcallers, callers, NULL);

    for (int i = 0; i < count; i++)
    {
        addr2line(callers[i]);
    }

    fprintf(stderr, "Call Stack:\n");
    for (int i = 0; i < count; i++)
    {
        printf("frame: %2d | 0x%" PRIX64 "\n", i, (uint64_t)callers[i]);
    }

    fprintf(stderr, "\n[DEBUG END]\n");
#endif

    char errmsg[ERRMSG_MAX];
    va_list ap;
    va_start(ap, format);
    vsnprintf(errmsg, ERRMSG_MAX, format, ap);
    va_end(ap);
    throw std::runtime_error(errmsg);
}

#define error(...)                                   \
    do                                               \
    {                                                \
        error_stop(__FILE__, __LINE__, __VA_ARGS__); \
    } while (0)

#endif

// DEBUG error message
#ifdef error
#undef error
#endif

#define error(...)                                   \
    do                                               \
    {                                                \
        error_stop(__FILE__, __LINE__, __VA_ARGS__); \
    } while (0)

#ifdef assert
#undef assert
#endif

/** redefine assert to use error and make backtraces */
#ifdef assert
#undef assert
#endif
#define assert(x)                              \
    do                                         \
    {                                          \
        if (!(x))                              \
            error("Assertion failed: %s", #x); \
    } while (0)

#ifdef _DEBUG

// this makes the following classes inspectable using lldb on OSX
// http://stackoverflow.com/questions/39680320/printing-debugging-libc-stl-with-xcode-lldb

#include <vector>
#include <list>
#include <string>

template class std::list<int>;
template class std::list<float>;
template class std::list<std::string>;

template class std::vector<int>;
template class std::vector<float>;
template class std::vector<std::string>;

#endif
