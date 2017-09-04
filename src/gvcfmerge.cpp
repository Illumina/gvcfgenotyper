#include "GVCFReader.hpp"

using namespace std;

int main(int argc, char **argv) 
{
    assert(argc>=2);
    string gvcf = argv[1];
    string ref = argv[2];
    GVCFReader g(gvcf,ref);
    return(0);
}
