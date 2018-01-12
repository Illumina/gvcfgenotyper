#include "GVCFMerger.hh"
#include <getopt.h>

static void usage()
{
    std::cerr << "\nAbout:   GVCF merging and genotyping for Illumina GVCFs" << std::endl;
    std::cerr << "Version: " << GIT_VERSION << std::endl;
    std::cerr << "Usage:   gvcfmerge -f ref.fa -l gvcf_list.txt" << std::endl;
    std::cerr << "" << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << "    -l, --list        <file>          plain text list of gvcfs to merge" << std::endl;
    std::cerr << "    -f, --fasta-ref   <file>          reference sequence" << std::endl;
    std::cerr << "    -o, --output-file <file>          output file name [stdout]" << std::endl;
    std::cerr
            << "    -O, --output-type <b|u|z|v>       b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]"
            << std::endl;
    std::cerr << "    -r, --regions     <region>        region to genotype eg. chr1 or chr20:5000000-6000000"
              << std::endl;
//    std::cerr << "    -@, --thread      INT             number of threads [0]" << std::endl;
    std::cerr << std::endl;
}

int main(int argc, char **argv)
{
    if (argc < 2)
    { usage(); }
    int c;
    string region = "";
    int n_threads = 0;
    string output_file = "";
    string output_type = "v";
    string gvcf_list = "";
    string reference_genome = "";
    
    static struct option loptions[] = {
            {"list",        1, 0, 'l'},
            {"fasta-ref",   1, 0, 'f'},
            {"output-file", 1, 0, 'o'},
            {"output-type", 1, 0, 'O'},
            {"region",      1, 0, 'r'},
            {"thread",      1, 0, '@'},

            {0,             0, 0, 0}
    };

    while ((c = getopt_long(argc, argv, "l:f:o:O:r:@:", loptions, NULL)) >= 0)
    {
        switch (c)
        {
            case 'l':
                gvcf_list = optarg;
                break;
            case 'f':
                reference_genome = optarg;
                break;
            case 'o':
                output_file = optarg;
                break;
            case 'O':
                output_type = optarg;
                break;
            case 'r':
                region = optarg;
                break;
            case '@':
                n_threads = atoi(optarg);
                break;
            default:
                if (optarg != NULL)
                { ggutils::die("Unknown argument:" + (string) optarg + "\n"); }
                else
                { ggutils::die("unrecognised argument"); }
        }
    }

    if (gvcf_list.empty())
    {
        ggutils::die("--list is required");
    }

    if (reference_genome.empty())
    {
        ggutils::die("--fasta-ref is required");
    }
    if (output_type != "b" && output_type != "z" && output_type != "v" && output_type != "u")
    {
        ggutils::die("invalid output type: " + output_type);
    }
    if (n_threads != 0)
    {
        ggutils::die("-@ is not implemented");
    }

    int buffer_size = 5000;
    std::vector<std::string> input_files;
    ggutils::read_text_file(gvcf_list, input_files);
    int is_file = 0;
    GVCFMerger g(input_files, output_file, output_type, reference_genome, buffer_size, region, is_file);

    g.write_vcf();
    return (0);
}
