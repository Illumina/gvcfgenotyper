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
#include <htslib/vcfutils.h>
#include "htslib/synced_bcf_reader.h"
}

int *zeros(int n);

typedef unsigned char byte;
typedef unsigned int uint;

bool fileexists(const string& fname);

template <typename T> void newMatrix(int nrow, int ncol, T** & out) { 
    out = new T * [nrow];
    for(int i = 0; i < nrow; i++) 
        out[i] = new T[ncol];  
    //  cout << out[nrow-1][ncol-1]<<endl;
};


template <typename T> void printMatrix(int nrow, int ncol, T* & out) { 
    for(int i = 0; i < nrow; i++) {
        for(int j = 0; j < ncol; j++) {
            cout << out[i * ncol + j] << " ";
        }
        cout << endl;
    }
    cout << endl;
};

template <typename T> void delMatrix(int nrow, int ncol, T** & out) { 
    for(int i = 0; i < nrow; i++) {
        delete[] out[i];
    }
    delete[] out;
};


//retruns first location where find_this occurs in_that (-1 if not found);
template <typename T> int find(T find_this, vector<T> in_that) { 
    for(int i = 0; i < in_that.size(); i++) {
        if(in_that[i] == find_this) {
            return(i);
        }
    }
    return(-1);
};


inline float argmax(float *P, int K, uint & maxi, uint & maxj) {
    //  assert(!isnan(*ptr));
    float *ptr = P;
    float maxval = *ptr;
    maxi=0;
    maxj=0;
    for(int i = 0; i < K; i++) {
        for(int j = 0; j < K; j++) {
            if(*ptr > maxval) {
                maxval = *ptr;
                maxi = i;
                maxj = j;
            }
            ptr++;
        }
    }
    return(maxval);
}

string percent(int num, int den);

inline void die(const string& s) {
    cerr << "ERROR: " << s << "\nExiting..." << endl;
    exit(1);
}

int strsplit(const string& input,const char split, vector<string> & out); 

string join(const vector<string>& input, const string& delim);

vector<int> match(const vector<string>& x, const vector<string>& y);

inline void warn(const string& s) {
    cerr << "WARNING: "<<s<<endl;
}


htsFile *vcf_wopen(const string& out_filename, string output_type);


int copyContigs(bcf_hdr_t *src,bcf_hdr_t *dst);


//simple dumps text from fname into output
int readTextFile(char *fname,vector<string> & output);
