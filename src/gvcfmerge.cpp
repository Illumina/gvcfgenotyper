#include "GVCFMerger.hpp"

using namespace std;

int main(int argc,char **argv)
{
    assert(argc>=2);
    vector<string> input_files;
    for(int i=2;i<argc;i++)
    {
    	input_files.push_back(argv[i]);
    }

    GVCFMerger g(input_files,"","v",argv[1],200);
    
    g.write_vcf();
    return(0);
}
