/**
 * @file metabolite_test.cpp
 * @brief Unit tests for the Metabolite class using Google Test
 */

#include <gtest/gtest.h>
#include <sstream>
#include "mosh/metabolite.h"

using namespace mosh;

/**
 * @class MetaboliteTest
 * @brief Test fixture for Metabolite class tests
 */
class MetaboliteTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create test metabolites
        glucose = new Metabolite("glc_D", "D-Glucose");
        atp = new Metabolite("atp", "Adenosine triphosphate");
    }

    void TearDown() override {
        delete glucose;
        delete atp;
    }

    Metabolite* glucose;
    Metabolite* atp;
};

// Test constructor and basic accessors
TEST_F(MetaboliteTest, ConstructorInitializesCorrectly) {
    EXPECT_EQ(glucose->name(), "glc_D");
    EXPECT_EQ(glucose->full_name(), "D-Glucose");
    EXPECT_FALSE(glucose->is_source());
    EXPECT_FALSE(glucose->is_cycle_met());
    EXPECT_FALSE(glucose->is_c_source());
    EXPECT_FALSE(glucose->is_dummy());
}

// Test name setters and getters
TEST_F(MetaboliteTest, NameSetterGetter) {
    glucose->set_name("glucose");
    EXPECT_EQ(glucose->name(), "glucose");
}

// Test index setter and getter
TEST_F(MetaboliteTest, IndexSetterGetter) {
    glucose->set_index(42);
    EXPECT_EQ(glucose->index(), 42);
}

// Test source flag
TEST_F(MetaboliteTest, SourceFlagSetterGetter) {
    EXPECT_FALSE(glucose->is_source());
    glucose->set_is_source(true);
    EXPECT_TRUE(glucose->is_source());
    glucose->set_is_source(false);
    EXPECT_FALSE(glucose->is_source());
}

// Test cycle metabolite flag
TEST_F(MetaboliteTest, CycleMetFlagSetterGetter) {
    EXPECT_FALSE(glucose->is_cycle_met());
    glucose->set_cycle_met(true);
    EXPECT_TRUE(glucose->is_cycle_met());
    glucose->set_cycle_met(false);
    EXPECT_FALSE(glucose->is_cycle_met());
}

// Test carbon source flag
TEST_F(MetaboliteTest, CarbonSourceFlagSetterGetter) {
    EXPECT_FALSE(glucose->is_c_source());
    glucose->set_c_source(true);
    EXPECT_TRUE(glucose->is_c_source());
    glucose->set_c_source(false);
    EXPECT_FALSE(glucose->is_c_source());
}

// Test dummy flag
TEST_F(MetaboliteTest, DummyFlagSetterGetter) {
    EXPECT_FALSE(glucose->is_dummy());
    glucose->set_is_dummy(true);
    EXPECT_TRUE(glucose->is_dummy());
    glucose->set_is_dummy(false);
    EXPECT_FALSE(glucose->is_dummy());
}

// Test write_to output
TEST_F(MetaboliteTest, WriteToOutputStream) {
    std::ostringstream oss;
    glucose->write_to(oss);
    std::string expected = "MET glc_D \"D-Glucose\"\n";
    EXPECT_EQ(oss.str(), expected);
}

// Test multiple flag combinations
TEST_F(MetaboliteTest, MultipleFlagCombinations) {
    glucose->set_is_source(true);
    glucose->set_c_source(true);
    glucose->set_cycle_met(true);

    EXPECT_TRUE(glucose->is_source());
    EXPECT_TRUE(glucose->is_c_source());
    EXPECT_TRUE(glucose->is_cycle_met());
    EXPECT_FALSE(glucose->is_dummy());
}

// Test independent metabolite instances
TEST_F(MetaboliteTest, IndependentInstances) {
    glucose->set_is_source(true);
    atp->set_is_source(false);

    EXPECT_TRUE(glucose->is_source());
    EXPECT_FALSE(atp->is_source());
    EXPECT_NE(glucose->name(), atp->name());
}

// Test empty name handling
TEST_F(MetaboliteTest, EmptyNameHandling) {
    Metabolite empty("", "");
    EXPECT_EQ(empty.name(), "");
    EXPECT_EQ(empty.full_name(), "");
}

// Test long name handling
TEST_F(MetaboliteTest, LongNameHandling) {
    std::string long_name(1000, 'x');
    Metabolite long_met(long_name, long_name);
    EXPECT_EQ(long_met.name(), long_name);
    EXPECT_EQ(long_met.full_name(), long_name);
}

// Test special characters in names
TEST_F(MetaboliteTest, SpecialCharactersInNames) {
    Metabolite special("h2o_c", "H2O[cytosol]");
    EXPECT_EQ(special.name(), "h2o_c");
    EXPECT_EQ(special.full_name(), "H2O[cytosol]");
}
