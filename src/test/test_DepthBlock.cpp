//
// Created by O'Connell, Jared on 10/16/17.
//

#include "test_helpers.hh"
#include "DepthBlock.hh"

TEST(DepthBlock, DepthBlock_def_constr)
{
    DepthBlock db;
    // nt _rid, _start, _end, _dp, _gq, _dpf;
    ASSERT_EQ(db._rid,0);
    ASSERT_EQ(db._start,0);
    ASSERT_EQ(db._end,0);
    ASSERT_EQ(db._dp,bcf_int32_missing);
    ASSERT_EQ(db._gq,bcf_int32_missing);
    ASSERT_EQ(db._dpf,bcf_int32_missing);
}

TEST(DepthBlock, DepthBlock_constr)
{
    DepthBlock db(5,15,20,10,11,12);
    ASSERT_EQ(db._rid,5);
    ASSERT_EQ(db._start,15);
    ASSERT_EQ(db._end,20);
    ASSERT_EQ(db._dp,10);
    ASSERT_EQ(db._dpf,11);
    ASSERT_EQ(db._gq,12);
}

TEST(DepthBlock, DepthBlock_intersect1)
{
    DepthBlock db1(5,1,100,10,10,10);
    DepthBlock db2(5,50,60,5,5,5);
    DepthBlock inters(5,50,60,10,10,10); 
    ASSERT_EQ(db1.intersect(db2),inters);
}

TEST(DepthBlock, DepthBlock_intersect2)
{
    DepthBlock db1(5,1,100,10,10,10);
    int rid=5, start=40, end=80;
    DepthBlock inters(5,40,80,10,10,10); 
    ASSERT_EQ(db1.intersect(rid,start,end),inters);
}

TEST(DepthBlock, DepthBlock_intersect_size1)
{
    DepthBlock db1(5,1,100,10,10,10);
    DepthBlock db2(5,50,60,5,5,5);
    ASSERT_EQ(db1.intersect_size(db2),11);

    DepthBlock db3(10,50,60,5,5,5);
    ASSERT_EQ(db1.intersect_size(db3),0);
}


TEST(DepthBlock, DepthBlock_intersect_size2)
{
    DepthBlock db1(5,1,100,10,10,10);
    int rid=5, start=40, end=80;
    ASSERT_EQ(db1.intersect_size(rid,start,end),41);
    rid = 11;
    ASSERT_EQ(db1.intersect_size(rid,start,end),0);
}

TEST(DepthBlock, DepthBlock_size)
{
    DepthBlock db1(5,1,100,10,10,10);
    ASSERT_EQ(db1.size(),100);
}

TEST(DepthBlock, DepthBlock_set_missing)
{
    DepthBlock db(5,1,100,10,10,10);
    db.set_missing();
    ASSERT_EQ(db._dp,bcf_int32_missing);
    ASSERT_EQ(db._gq,bcf_int32_missing);
    ASSERT_EQ(db._dpf,bcf_int32_missing);
}

TEST(DepthBlock, DepthBlock_set_zero)
{
    DepthBlock db(5,1,100,10,10,10);
    db.zero();
    ASSERT_EQ(db._dp,0);
    ASSERT_EQ(db._gq,0);
    ASSERT_EQ(db._dpf,0);
}


TEST(DepthBlock, DepthBlock_add1)
{
    DepthBlock db1(5,1,100,10,10,10);
    DepthBlock db2(5,101,200,20,20,20);
    db1.add(db2);
    ASSERT_EQ(db1._rid,5);
    ASSERT_EQ(db1._start,1);
    ASSERT_EQ(db1._end,200);
    ASSERT_EQ(db1._dp,15);
    ASSERT_EQ(db1._dpf,15);
    ASSERT_EQ(db1._gq,15);
}

TEST(DepthBlock, DepthBlock_add2)
{
    DepthBlock db1(5,1,100,10,10,10);
    DepthBlock db2(5,101,150,20,20,20);
    db1.add(db2);
    ASSERT_EQ(db1._rid,5);
    ASSERT_EQ(db1._start,1);
    ASSERT_EQ(db1._end,150);
    ASSERT_EQ(db1._dp,13);
    ASSERT_EQ(db1._dpf,13);
    ASSERT_EQ(db1._gq,13);
}
