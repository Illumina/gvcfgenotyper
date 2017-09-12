#include "GVCFMerger.hpp"

using namespace std;

int main(int argc,char **argv)
{
    if(argc<2)
    {
	cerr << "Usage: gvcfmerge reference.fa s1.genome.vcf.gz s2.genome.vcf.gz ... sN.genome.vcf.gz" <<endl;
	return(0);
    }
    vector<string> input_files;
    for(int i=2;i<argc;i++)
    {
    	input_files.push_back(argv[i]);
    }
    int buffer_size = 200;
    GVCFMerger g(input_files,"","v",argv[1],buffer_size);
    
    g.write_vcf();
    return(0);
}
