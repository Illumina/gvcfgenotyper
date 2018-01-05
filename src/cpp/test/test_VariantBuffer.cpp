//
// Created by Schulz-Trieglaff, Ole on 11/03/17.
//

#include "test_helpers.hh"
#include "VariantBuffer.hh"

TEST(VariantBuffer, VariantBuffer_def_constr)
{
    VariantBuffer vb;
    ASSERT_EQ(vb.Size(),(size_t)0);
    ASSERT_EQ(vb.IsEmpty(),true);
}

TEST(VariantBuffer, VariantBuffer_push_back)
{
    VariantBuffer vb;

    auto hdr = get_header();
    int rid=13;
    int pos=600;
    auto rec1 = generate_record(hdr,rid,pos,"A,T");
    vb.PushBack(hdr, rec1);
    ASSERT_EQ(vb.Size(),(size_t)1);
    ASSERT_EQ(vb.IsEmpty(),false);

    rid=13;
    pos=800;
    auto rec2 = generate_record(hdr,rid,pos,"A,G");
    vb.PushBack(hdr, rec2);
    ASSERT_EQ(vb.Size(),(size_t)2);
    ASSERT_EQ(vb.IsEmpty(),false);
}

TEST(VariantBuffer, VariantBuffer_flush_buffer1)
{
    VariantBuffer vb;

    auto hdr = get_header();
    auto rec1 = generate_record(hdr,13,600,"A,T");
    vb.PushBack(hdr, rec1);
    auto rec2 = generate_record(hdr,13,800,"A,G");
    vb.PushBack(hdr, rec2);
    auto rec3 = generate_record(hdr,13,1000,"A,TTG");
    vb.PushBack(hdr, rec3);

    ASSERT_EQ(vb.FlushBuffer(13, 800),2);
    ASSERT_EQ(vb.Size(),(size_t)1);
    ASSERT_EQ(vb.IsEmpty(),false);

    ASSERT_EQ(vb.FlushBuffer(13, 1500),1);
    ASSERT_EQ(vb.Size(),(size_t)0);
    ASSERT_EQ(vb.IsEmpty(),true);
}

TEST(VariantBuffer, VariantBuffer_flush_buffer2)
{
    VariantBuffer vb;

    auto hdr = get_header();
    auto rec1 = generate_record(hdr,13,600,"A,T");
    vb.PushBack(hdr, rec1);
    auto rec2 = generate_record(hdr,13,800,"A,G");
    vb.PushBack(hdr, rec2);
    auto rec3 = generate_record(hdr,13,1000,"A,T");
    vb.PushBack(hdr, rec3);

    auto rec2a = generate_record(hdr,13,800,"A,G");
    ASSERT_EQ(vb.FlushBuffer(rec2a),2);
    ASSERT_EQ(vb.Size()==1,true);
    ASSERT_EQ(vb.IsEmpty(),false);

    auto rec3a = generate_record(hdr,13,1000,"A,T");
    ASSERT_EQ(vb.FlushBuffer(rec3a),1);
    ASSERT_EQ(vb.Size()==0,true);
    ASSERT_EQ(vb.IsEmpty(),true);
}

TEST(VariantBuffer, VariantBuffer_flush_buffer3)
{
    VariantBuffer vb;

    auto hdr = get_header();
    auto rec1 = generate_record(hdr,13,600,"A,T");
    vb.PushBack(hdr, rec1);
    auto rec2 = generate_record(hdr,13,800,"A,G");
    vb.PushBack(hdr, rec2);
    auto rec3 = generate_record(hdr,13,1000,"A,TTG");
    vb.PushBack(hdr, rec3);

    ASSERT_EQ(vb.FlushBuffer(),3);
    ASSERT_EQ(vb.Size()==0,true);
    ASSERT_EQ(vb.IsEmpty(),true);
}

TEST(VariantBuffer, VariantBuffer_size)
{
    VariantBuffer vb;
    ASSERT_EQ(vb.Size()==0,true);
    auto hdr = get_header();
    auto rec1 = generate_record(hdr,13,600,"A,T");
    vb.PushBack(hdr, rec1);
    ASSERT_EQ(vb.Size()==1,true);
    auto rec2 = generate_record(hdr,13,800,"A,G");
    vb.PushBack(hdr, rec2);
    ASSERT_EQ(vb.Size()==2,true);
}

TEST(VariantBuffer, VariantBuffer_empty)
{
    VariantBuffer vb;
    ASSERT_EQ(vb.IsEmpty(),true);
    auto hdr = get_header();
    auto rec1 = generate_record(hdr,13,600,"A,T");
    vb.PushBack(hdr, rec1);
    ASSERT_EQ(vb.IsEmpty(),false);
}

TEST(VariantBuffer, VariantBuffer_back)
{
    VariantBuffer vb;
    auto hdr = get_header();
    auto rec1 = generate_record(hdr,13,600,"A,T");
    vb.PushBack(hdr, rec1);
    ASSERT_EQ(vb.Back(),rec1);
}

TEST(VariantBuffer, VariantBuffer_push_back_duplicate)
{
    VariantBuffer vb;
    auto hdr = get_header();
    auto rec1 = generate_record(hdr,20,20000,"A,T");
    vb.PushBack(hdr, rec1);
    auto rec2 = generate_record(hdr,20,20000,"A,T");
    vb.PushBack(hdr, rec2);
    ASSERT_EQ(vb.Size(),(size_t)1);
    ASSERT_EQ(vb.Back(),rec1);
    ASSERT_EQ(vb.GetNumDuplicatedRecords(),(size_t)1);
    vb.FlushBuffer();
}

TEST(VariantBuffer, VariantBuffer_push_back_duplicate_hom_ref)
{
    VariantBuffer vb;
    auto hdr = get_header();
    bcf1_t* rec1 = generate_record(hdr,20,20000,"A,T","0/0");
    vb.PushBack(hdr, rec1);
    bcf1_t* rec2 = generate_record(hdr,20,20000,"A,T","0/1");
    // need a copy of rec2 here since rec2 will be destroyed
    bcf1_t* rec2b = bcf_dup(rec2);
    vb.PushBack(hdr, rec2);

    ASSERT_EQ(vb.Size(),(size_t)1);
    bcf1_t* rec3 = vb.Pop();
    bool cmp = ggutils::bcf1_equal(rec2b,rec3) && !ggutils::is_hom_ref(hdr,rec3);
    ASSERT_EQ(cmp,true);

    ASSERT_EQ(vb.GetNumDuplicatedRecords(),(size_t)1);
    vb.FlushBuffer();
}
