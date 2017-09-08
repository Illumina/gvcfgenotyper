#include "GVCFReader.hpp"

int main(int argc,char **argv)
{
    cerr << "GVCF: "<<argv[1]<<endl;
    cerr << "Reference: "<<argv[2]<<endl;
    assert(argc>=2);
    string gvcf = argv[1];
    string ref = argv[2];
    GVCFReader g(gvcf,ref,200);
    int num_lines_to_read = 100;

//    DepthBlock db(0,0,0,0,0,0);
    cerr << integer_thing(333)<<endl;
    const    bcf_hdr_t *hdr = g.getHeader();    
    bcf1_t *line = g.pop();
    while(line!=NULL)
    {
	cout << bcf_hdr_id2name(hdr,line->rid)<<":"<<line->pos+1<<":"<<line->d.allele[0]<<":"<<line->d.allele[1]<<endl;
	assert(line->n_allele==2);
	bcf_destroy(line);
	line = g.pop();
    }
    assert(g.empty());
    return(0);
}
