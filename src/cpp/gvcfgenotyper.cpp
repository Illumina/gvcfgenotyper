#include "GVCFMerger.hh"
#include <getopt.h>

#include "spdlog.h"

static void usage()
{
    std::cerr << "\nAbout:   GVCF merging and genotyping for Illumina GVCFs" << std::endl;
    std::cerr << "Version: " << GG_VERSION << std::endl;
    std::cerr << "Usage:   gvcfmerge -f ref.fa -l gvcf_list.txt" << std::endl;
    std::cerr << "" << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << "    -l, --list          <file>          plain text list of gvcfs to merge" << std::endl;
    std::cerr << "    -f, --fasta-ref     <file>          reference sequence" << std::endl;
    std::cerr << "    -o, --output-file   <file>          output file name [stdout]" << std::endl;
    std::cerr
            << "    -O, --output-type     <b|u|z|v>       b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]"
            << std::endl;
    std::cerr << "    -r, --region        <region>        region to genotype eg. chr1 or chr20:5000000-6000000"
              << std::endl;
    std::cerr << "    -M, --max-alleles   <file>          maximum number of alleles [20]" << std::endl;
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
    string log_file = "gvcfgenotyper."+ggutils::string_time()+"."+to_string(getpid())+".log";
    string output_type = "v";
    string gvcf_list = "";
    string reference_genome = "";
    size_t max_alleles = 20;

    //This is a hidden flag that when true will drop variants with reference mismatches rather than exit with error (this is ill advised).
    bool ignore_non_matching_ref=false;
    // Another hidden flag to force processing of gvcf files with duplicate sample names
    bool force_samples=false;

    static struct option loptions[] = {
            {"list",        1, 0, 'l'},
            {"fasta-ref",   1, 0, 'f'},
            {"output-file", 1, 0, 'o'},
            {"output-type", 1, 0, 'O'},
            {"log-file",    1, 0, 'g'},
            {"region",      1, 0, 'r'},
            {"thread",      1, 0, '@'},
            {"max-alleles", 1, 0, 'M'},
	        {"ignore-non-matching-ref",0,0,1},
	        {"force-samples",0,0,'s'},
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
            case 'M':
                output_type = optarg;
                break;
            case 'g':
                log_file = optarg;
                break;
            case 'r':
                region = optarg;
                break;
            case '@':
                n_threads = atoi(optarg);
                break;
            case 1:
	            ignore_non_matching_ref=true;break;
            case 's':
	            force_samples=true;
                break;
	        default:
	            if (optarg != NULL)
		            ggutils::die("Unknown argument:" + (string) optarg + "\n");
	            else
		            ggutils::die("unrecognised argument");
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
    std::cerr << "Logging output to " <<log_file<<std::endl;

    // register logger, name of outfile can be set by user on the cmd line
    std::shared_ptr<spdlog::logger> lg = spdlog::basic_logger_mt("gg_logger", log_file);
    // format: "*** [YYYY-MM-DD HH:MM:SS] [thread] [loglevel] message ***"
    spdlog::set_pattern(" [%c] [%l] %v");
    std::string commandline = argv[0];
    for(int i=1;i<argc;i++) commandline += (" " + (string)argv[i]);
    lg->info("Command line: "+commandline);
    lg->info("Starting GVCF merging");

    int buffer_size = 5000;
    std::vector<std::string> input_files;
    ggutils::read_text_file(gvcf_list, input_files);
    int is_file = 0;
    GVCFMerger g(input_files, output_file, output_type, reference_genome, buffer_size, region, is_file, ignore_non_matching_ref, force_samples);
    g.SetMaxAlleles(max_alleles);
    g.write_vcf();

    lg->info("Done");
    spdlog::drop_all();
    return (0);
}
