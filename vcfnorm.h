#include "bcftools.h"
#include "rbuf.h"
#include <htslib/faidx.h>


// for -m+, mapping from allele indexes of a single input record
// to allele indexes of output record
typedef struct
{
    int nals, mals, *map;
}
map_t;

typedef struct
{
    char *tseq, *seq;
    int mseq;
    bcf1_t **lines, **tmp_lines, **alines, **blines, *mrow_out;
    int ntmp_lines, mtmp_lines, nalines, malines, nblines, mblines;
    map_t *maps;     // mrow map for each buffered record
    char **als;
    int mmaps, nals, mals;
    uint8_t *tmp_arr1, *tmp_arr2, *diploid;
    int ntmp_arr1, ntmp_arr2;
    kstring_t *tmp_str;
    kstring_t *tmp_als, tmp_als_str;
    int ntmp_als;
    rbuf_t rbuf;
    int buf_win;            // maximum distance between two records to consider
    int aln_win;            // the realignment window size (maximum repeat size)
    bcf_srs_t *files;       // using the synced reader only for -r option
    bcf_hdr_t *hdr;
    faidx_t *fai;
    struct { int tot, set, swap; } nref;
    char **argv, *output_fname, *ref_fname, *vcf_fname, *region, *targets;
    int argc, rmdup, output_type, n_threads, check_ref, strict_filter, do_indels;
    int nchanged, nskipped, nsplit, ntotal, mrows_op, mrows_collapse, parsimonious;
}
args_t;

args_t *init_vcfnorm(  bcf_hdr_t *hdr,char *ref);


void split_multiallelic_to_biallelics(args_t *args, bcf1_t *line);
void write_vcf_line(args_t *args,bcf1_t *line,htsFile *out);
//void flush_buffer(args_t *args, htsFile *file, int n);
void destroy_data(args_t *args);
int realign(args_t *args, bcf1_t *line);
