#include "GVCFReader.cpp"

using namespace std;

int main(int argc, char **argv) 
{
    assert(argc>=2);
    GVCFReader g(argv[1],argv[2]);
    return(0);
}
