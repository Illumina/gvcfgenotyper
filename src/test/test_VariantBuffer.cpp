//
// Created by Schulz-Trieglaff, Ole on 11/03/17.
//

#include "test_helpers.hh"
#include "VariantBuffer.hh"

TEST(VariantBuffer, VariantBuffer_def_constr)
{
    VariantBuffer vb;
    ASSERT_EQ(vb.size(),(size_t)0);
    ASSERT_EQ(vb.empty(),true);
}

TEST(VariantBuffer, VariantBuffer_push_back)
{
    VariantBuffer vb;

    auto hdr = get_header();
    int rid=13;
    int pos=600;
    auto rec1 = generate_record(hdr,rid,pos,"A,T");
    vb.push_back(hdr,rec1);
    ASSERT_EQ(vb.size(),(size_t)1);
    ASSERT_EQ(vb.empty(),false);

    rid=13;
    pos=800;
    auto rec2 = generate_record(hdr,rid,pos,"A,G");
    vb.push_back(hdr,rec2);
    ASSERT_EQ(vb.size(),(size_t)2);
    ASSERT_EQ(vb.empty(),false);
}

TEST(VariantBuffer, VariantBuffer_flush_buffer1)
{
    VariantBuffer vb;

    auto hdr = get_header();
    auto rec1 = generate_record(hdr,13,600,"A,T");
    vb.push_back(hdr,rec1);
    auto rec2 = generate_record(hdr,13,800,"A,G");
    vb.push_back(hdr,rec2);
    auto rec3 = generate_record(hdr,13,1000,"A,TTG");
    vb.push_back(hdr,rec3);

    ASSERT_EQ(vb.flush_buffer(13,800),2);
    ASSERT_EQ(vb.size(),(size_t)1);
    ASSERT_EQ(vb.empty(),false);

    ASSERT_EQ(vb.flush_buffer(13,1500),1);
    ASSERT_EQ(vb.size(),(size_t)0);
    ASSERT_EQ(vb.empty(),true);
}

TEST(VariantBuffer, VariantBuffer_flush_buffer2)
{
    VariantBuffer vb;

    auto hdr = get_header();
    auto rec1 = generate_record(hdr,13,600,"A,T");
    vb.push_back(hdr,rec1);
    auto rec2 = generate_record(hdr,13,800,"A,G");
    vb.push_back(hdr,rec2);
    auto rec3 = generate_record(hdr,13,1000,"A,T");
    vb.push_back(hdr,rec3);

    auto rec2a = generate_record(hdr,13,800,"A,G");
    ASSERT_EQ(vb.flush_buffer(rec2a),2);
    ASSERT_EQ(vb.size()==1,true);
    ASSERT_EQ(vb.empty(),false);

    auto rec3a = generate_record(hdr,13,1000,"A,T");
    ASSERT_EQ(vb.flush_buffer(rec3a),1);
    ASSERT_EQ(vb.size()==0,true);
    ASSERT_EQ(vb.empty(),true);
}

TEST(VariantBuffer, VariantBuffer_flush_buffer3)
{
    VariantBuffer vb;

    auto hdr = get_header();
    auto rec1 = generate_record(hdr,13,600,"A,T");
    vb.push_back(hdr,rec1);
    auto rec2 = generate_record(hdr,13,800,"A,G");
    vb.push_back(hdr,rec2);
    auto rec3 = generate_record(hdr,13,1000,"A,TTG");
    vb.push_back(hdr,rec3);

    ASSERT_EQ(vb.flush_buffer(),3);
    ASSERT_EQ(vb.size()==0,true);
    ASSERT_EQ(vb.empty(),true);
}

TEST(VariantBuffer, VariantBuffer_size)
{
    VariantBuffer vb;
    ASSERT_EQ(vb.size()==0,true);
    auto hdr = get_header();
    auto rec1 = generate_record(hdr,13,600,"A,T");
    vb.push_back(hdr,rec1);
    ASSERT_EQ(vb.size()==1,true);
    auto rec2 = generate_record(hdr,13,800,"A,G");
    vb.push_back(hdr,rec2);
    ASSERT_EQ(vb.size()==2,true);
}

TEST(VariantBuffer, VariantBuffer_empty)
{
    VariantBuffer vb;
    ASSERT_EQ(vb.empty(),true);
    auto hdr = get_header();
    auto rec1 = generate_record(hdr,13,600,"A,T");
    vb.push_back(hdr,rec1);
    ASSERT_EQ(vb.empty(),false);
}

TEST(VariantBuffer, VariantBuffer_back)
{
    VariantBuffer vb;
    auto hdr = get_header();
    auto rec1 = generate_record(hdr,13,600,"A,T");
    vb.push_back(hdr,rec1);
    ASSERT_EQ(vb.back(),rec1);
}

TEST(VariantBuffer, VariantBuffer_push_back_duplicate)
{
    VariantBuffer vb;
    auto hdr = get_header();
    auto rec1 = generate_record(hdr,20,20000,"A,T");
    vb.push_back(hdr,rec1);
    auto rec2 = generate_record(hdr,20,20000,"A,T");
    vb.push_back(hdr,rec2);
    ASSERT_EQ(vb.size(),(size_t)1);
    ASSERT_EQ(vb.back(),rec1);
    ASSERT_EQ(vb.get_num_duplicated_records(),(size_t)1);
    vb.flush_buffer();
}

TEST(VariantBuffer, VariantBuffer_push_back_duplicate_hom_ref)
{
    VariantBuffer vb;
    auto hdr = get_header();
    bcf1_t* rec1 = generate_record(hdr,20,20000,"A,T","0/0");
    vb.push_back(hdr,rec1);
    bcf1_t* rec2 = generate_record(hdr,20,20000,"A,T","0/1");
    // need a copy of rec2 here since rec2 will be destroyed
    bcf1_t* rec2b = bcf_dup(rec2);
    vb.push_back(hdr,rec2);

    ASSERT_EQ(vb.size(),(size_t)1);
    bcf1_t* rec3 = vb.pop();
    bool cmp = bcf1_equal(rec2b,rec3) && !is_hom_ref(hdr,rec3);
    ASSERT_EQ(cmp,true);

    ASSERT_EQ(vb.get_num_duplicated_records(),(size_t)1);
    vb.flush_buffer();
}
