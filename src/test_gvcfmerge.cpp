#include "gtest/gtest.h"
#include "common.hpp"
#include <libgen.h>

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

GTestEnvironment* g_testenv = nullptr;

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    // determine base path
    char* bn = dirname(argv[0]);
    char actualpath[PATH_MAX + 1];
    realpath(bn, actualpath);

    ::testing::AddGlobalTestEnvironment(g_testenv = new GTestEnvironment(actualpath));
    int res = RUN_ALL_TESTS();

    return res;
}
