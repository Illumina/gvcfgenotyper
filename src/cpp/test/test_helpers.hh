//
// Created by O'Connell, Jared on 10/23/17.
//

#ifndef TEST_HELPERS_H
#define TEST_HELPERS_H

extern "C" {
#include <htslib/hts.h>
#include <htslib/vcf.h>
}

#include <dirent.h>
#include "gtest/gtest.h"
#include "common.hh"
#include <string>


bcf1_t *generate_record(bcf_hdr_t *hdr,int rid,int pos,const std::string & alleles, const std::string& gt);

bcf1_t *generate_record(bcf_hdr_t *hdr,int rid,int pos,const std::string & alleles);

bcf1_t *generate_record(bcf_hdr_t *hdr,const std::string & vcfrow);

bcf1_t *generate_record(bcf_hdr_t *hdr,const kstring_t & vcfrow);

bcf_hdr_t *get_header();

void update_record(bcf_hdr_t *hdr,int rid,int pos,const std::string & alleles,bcf1_t *record);


#endif
