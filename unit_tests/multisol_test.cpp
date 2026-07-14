/**
 * @file multisol_test.cpp
 * @brief Unit tests for the MultiSol class using Google Test
 */

#include <gtest/gtest.h>
#include "mosh/multisol.h"
#include "mosh/solution.h"
#include "mosh/scenario.h"
#include "mosh/params.h"
#include "mosh/metabolite.h"
#include "mosh/reaction.h"

using namespace mosh;

/**
 * @class MultiSolTest
 * @brief Test fixture for MultiSol class tests
 */
class MultiSolTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create scenario
        scenario = new Scenario();

        // Add metabolites
        auto glucose = std::make_shared<Metabolite>("glc_D", "D-Glucose");
        auto atp = std::make_shared<Metabolite>("atp", "ATP");
        scenario->add_metabolite(glucose);
        scenario->add_metabolite(atp);

        // Add reactions
        auto rxn1 = std::make_shared<Reaction>("R1", 1.0, 1000.0, true, true);
        auto rxn2 = std::make_shared<Reaction>("R2", 2.0, 500.0, true, true);
        scenario->add_reaction(rxn1);
        scenario->add_reaction(rxn2);

        // Create params
        params = new Params();

        // Create MultiSol
        multisol = new MultiSol(scenario, params);
    }

    void TearDown() override {
        delete multisol;
        delete params;
        delete scenario;
    }

    Scenario* scenario;
    Params* params;
    MultiSol* multisol;
};

// Test constructor
TEST_F(MultiSolTest, ConstructorInitializesCorrectly) {
    EXPECT_EQ(multisol->scenario(), scenario);
    EXPECT_EQ(multisol->params(), params);
}

// Test initial solution count
TEST_F(MultiSolTest, InitialSolutionCount) {
    EXPECT_EQ(multisol->num_solutions(), 0);
}

// Test scenario accessor
TEST_F(MultiSolTest, ScenarioAccessor) {
    EXPECT_EQ(multisol->scenario(), scenario);
}

// Test params accessor
TEST_F(MultiSolTest, ParamsAccessor) {
    EXPECT_EQ(multisol->params(), params);
}

// Test empty MultiSol
TEST_F(MultiSolTest, EmptyMultiSol) {
    EXPECT_EQ(multisol->num_solutions(), 0);
}
