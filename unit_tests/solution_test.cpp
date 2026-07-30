/**
 * @file solution_test.cpp
 * @brief Unit tests for the Solution class using Google Test
 */

#include <gtest/gtest.h>
#include "mosh/solution.h"
#include "mosh/scenario.h"
#include "mosh/params.h"
#include "mosh/metabolite.h"
#include "mosh/reaction.h"

using namespace mosh;

/**
 * @class SolutionTest
 * @brief Test fixture for Solution class tests
 */
class SolutionTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create scenario
        scenario = new Scenario();

        // Add metabolites
        auto glucose = std::make_shared<Metabolite>("glc_D", "D-Glucose");
        auto atp = std::make_shared<Metabolite>("atp", "ATP");
        auto adp = std::make_shared<Metabolite>("adp", "ADP");
        scenario->add_metabolite(glucose);
        scenario->add_metabolite(atp);
        scenario->add_metabolite(adp);

        // Add reactions
        auto rxn1 = std::make_shared<Reaction>("R1", 1.0, 1000.0, true, true);
        auto rxn2 = std::make_shared<Reaction>("R2", 2.0, 500.0, true, true);
        auto rxn3 = std::make_shared<Reaction>("R3", 3.0, 250.0, true, true);
        scenario->add_reaction(rxn1);
        scenario->add_reaction(rxn2);
        scenario->add_reaction(rxn3);

        // Create params
        params = new Params();

        // Create solution
        solution = new Solution(scenario, params);
    }

    void TearDown() override {
        delete solution;
        delete params;
        delete scenario;
    }

    Scenario* scenario;
    Params* params;
    Solution* solution;
};

// Test constructor
TEST_F(SolutionTest, ConstructorInitializesCorrectly) {
    EXPECT_EQ(solution->scenario(), scenario);
    EXPECT_EQ(solution->params(), params);
}

// Test flux initialization to zero
TEST_F(SolutionTest, FluxInitializedToZero) {
    for (size_t i = 0; i < scenario->num_reactions(); ++i) {
        EXPECT_EQ(solution->flux(i), 0.0);
    }
}

// Test set and get flux by index
TEST_F(SolutionTest, SetGetFluxByIndex) {
    solution->set_flux(0, 5.5);
    solution->set_flux(1, 10.2);
    solution->set_flux(2, 3.7);

    EXPECT_EQ(solution->flux(0), 5.5);
    EXPECT_EQ(solution->flux(1), 10.2);
    EXPECT_EQ(solution->flux(2), 3.7);
}

// Test set and get flux by reaction pointer
TEST_F(SolutionTest, SetGetFluxByReaction) {
    const Reaction* rxn1 = scenario->reaction(0);
    const Reaction* rxn2 = scenario->reaction(1);

    solution->set_flux(rxn1, 7.3);
    solution->set_flux(rxn2, 4.8);

    EXPECT_EQ(solution->flux(rxn1), 7.3);
    EXPECT_EQ(solution->flux(rxn2), 4.8);
}

// Test uses_react with positive flux
TEST_F(SolutionTest, UsesReactPositiveFlux) {
    solution->set_flux(0, 5.0);
    EXPECT_TRUE(solution->uses_react(0));
}

// Test uses_react with zero flux
TEST_F(SolutionTest, UsesReactZeroFlux) {
    solution->set_flux(0, 0.0);
    EXPECT_FALSE(solution->uses_react(0));
}

// Test uses_react with negative flux (should be false)
TEST_F(SolutionTest, UsesReactNegativeFlux) {
    solution->set_flux(0, -5.0);
    EXPECT_FALSE(solution->uses_react(0));
}

// Test uses_react by reaction pointer
TEST_F(SolutionTest, UsesReactByReactionPointer) {
    const Reaction* rxn = scenario->reaction(0);
    solution->set_flux(rxn, 10.0);

    EXPECT_TRUE(solution->uses_react(rxn));
}

// Test fill function
TEST_F(SolutionTest, FillAllFluxes) {
    solution->fill(42.0);

    for (size_t i = 0; i < scenario->num_reactions(); ++i) {
        EXPECT_EQ(solution->flux(i), 42.0);
    }
}

// Test copy constructor
TEST_F(SolutionTest, CopyConstructor) {
    solution->set_flux(0, 5.0);
    solution->set_flux(1, 10.0);

    Solution copy(solution);

    EXPECT_EQ(copy.flux(0), 5.0);
    EXPECT_EQ(copy.flux(1), 10.0);
    EXPECT_EQ(copy.scenario(), scenario);
    EXPECT_EQ(copy.params(), params);
}

// Test multiple flux updates
TEST_F(SolutionTest, MultipleFluxUpdates) {
    solution->set_flux(0, 1.0);
    EXPECT_EQ(solution->flux(0), 1.0);

    solution->set_flux(0, 2.0);
    EXPECT_EQ(solution->flux(0), 2.0);

    solution->set_flux(0, 0.0);
    EXPECT_EQ(solution->flux(0), 0.0);
}

// Test flux values are independent
TEST_F(SolutionTest, FluxValuesIndependent) {
    solution->set_flux(0, 5.0);
    solution->set_flux(1, 10.0);
    solution->set_flux(2, 15.0);

    EXPECT_EQ(solution->flux(0), 5.0);
    EXPECT_EQ(solution->flux(1), 10.0);
    EXPECT_EQ(solution->flux(2), 15.0);
}

// Test large flux values
TEST_F(SolutionTest, LargeFluxValues) {
    solution->set_flux(0, 1e6);
    EXPECT_EQ(solution->flux(0), 1e6);

    solution->set_flux(1, 1e-6);
    EXPECT_EQ(solution->flux(1), 1e-6);
}

// Test negative flux values
TEST_F(SolutionTest, NegativeFluxValues) {
    solution->set_flux(0, -100.0);
    EXPECT_EQ(solution->flux(0), -100.0);
}

// Test scenario accessor
TEST_F(SolutionTest, ScenarioAccessor) {
    EXPECT_EQ(solution->scenario(), scenario);
    EXPECT_EQ(solution->scenario()->num_reactions(), 3);
    EXPECT_EQ(solution->scenario()->num_metabolites(), 3);
}

// Test params accessor
TEST_F(SolutionTest, ParamsAccessor) {
    EXPECT_EQ(solution->params(), params);
}
