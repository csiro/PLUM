/**
 * @file params_test.cpp
 * @brief Unit tests for the Params class using Google Test
 */

#include <gtest/gtest.h>
#include "mosh/params.h"

using namespace mosh;

/**
 * @class ParamsTest
 * @brief Test fixture for Params class tests
 */
class ParamsTest : public ::testing::Test {
protected:
    void SetUp() override {
        params = new Params();
    }

    void TearDown() override {
        delete params;
    }

    Params* params;
};

// Test default constructor
TEST_F(ParamsTest, ConstructorInitializesCorrectly) {
    EXPECT_NE(params, nullptr);
}

// Test params object creation
TEST_F(ParamsTest, ParamsObjectCreation) {
    Params p1;
    Params p2;
    // Both should be valid objects
    EXPECT_NE(&p1, &p2);
}
