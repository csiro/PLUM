/**
 * @file scenario_test.cpp
 * @brief Unit tests for the Scenario class using Google Test
 */

#include <gtest/gtest.h>
#include "mosh/scenario.h"
#include "mosh/metabolite.h"
#include "mosh/reaction.h"

using namespace mosh;

/**
 * @class ScenarioTest
 * @brief Test fixture for Scenario class tests
 */
class ScenarioTest : public ::testing::Test {
protected:
    void SetUp() override {
        scenario = new Scenario();
    }

    void TearDown() override {
        delete scenario;
    }

    Scenario* scenario;
};

// Test default constructor
TEST_F(ScenarioTest, ConstructorInitializesCorrectly) {
    EXPECT_EQ(scenario->num_metabolites(), 0);
    EXPECT_EQ(scenario->num_reactions(), 0);
}

// Test adding metabolites
TEST_F(ScenarioTest, AddMetabolite) {
    auto glucose = std::make_shared<Metabolite>("glc_D", "D-Glucose");
    scenario->add_metabolite(glucose);

    EXPECT_EQ(scenario->num_metabolites(), 1);
    EXPECT_EQ(scenario->metabolite(0)->name(), "glc_D");
    EXPECT_EQ(scenario->metabolite(0)->index(), 0);
}

// Test adding multiple metabolites
TEST_F(ScenarioTest, AddMultipleMetabolites) {
    auto glucose = std::make_shared<Metabolite>("glc_D", "D-Glucose");
    auto atp = std::make_shared<Metabolite>("atp", "ATP");
    auto adp = std::make_shared<Metabolite>("adp", "ADP");

    scenario->add_metabolite(glucose);
    scenario->add_metabolite(atp);
    scenario->add_metabolite(adp);

    EXPECT_EQ(scenario->num_metabolites(), 3);
    EXPECT_EQ(scenario->metabolite(0)->name(), "glc_D");
    EXPECT_EQ(scenario->metabolite(1)->name(), "atp");
    EXPECT_EQ(scenario->metabolite(2)->name(), "adp");
}

// Test metabolite index assignment
TEST_F(ScenarioTest, MetaboliteIndexAssignment) {
    auto met1 = std::make_shared<Metabolite>("met1", "Metabolite 1");
    auto met2 = std::make_shared<Metabolite>("met2", "Metabolite 2");

    scenario->add_metabolite(met1);
    scenario->add_metabolite(met2);

    EXPECT_EQ(met1->index(), 0);
    EXPECT_EQ(met2->index(), 1);
}

// Test adding reactions
TEST_F(ScenarioTest, AddReaction) {
    auto rxn = std::make_shared<Reaction>("R1", 1.0, 1000.0, true, true);
    scenario->add_reaction(rxn);

    EXPECT_EQ(scenario->num_reactions(), 1);
    EXPECT_EQ(scenario->reaction(0)->name(), "R1");
}

// Test adding multiple reactions
TEST_F(ScenarioTest, AddMultipleReactions) {
    auto rxn1 = std::make_shared<Reaction>("R1", 1.0, 1000.0, true, true);
    auto rxn2 = std::make_shared<Reaction>("R2", 2.0, 500.0, true, true);
    auto rxn3 = std::make_shared<Reaction>("R3", 3.0, 250.0, true, true);

    scenario->add_reaction(rxn1);
    scenario->add_reaction(rxn2);
    scenario->add_reaction(rxn3);

    EXPECT_EQ(scenario->num_reactions(), 3);
    EXPECT_EQ(scenario->reaction(0)->name(), "R1");
    EXPECT_EQ(scenario->reaction(1)->name(), "R2");
    EXPECT_EQ(scenario->reaction(2)->name(), "R3");
}

// Test reaction index assignment
TEST_F(ScenarioTest, ReactionIndexAssignment) {
    auto rxn1 = std::make_shared<Reaction>("R1", 1.0, 1000.0, true, true);
    auto rxn2 = std::make_shared<Reaction>("R2", 2.0, 500.0, true, true);

    scenario->add_reaction(rxn1);
    scenario->add_reaction(rxn2);

    EXPECT_EQ(rxn1->index(), 0);
    EXPECT_EQ(rxn2->index(), 1);
}

// Test metabolites accessor
TEST_F(ScenarioTest, MetabolitesAccessor) {
    auto met1 = std::make_shared<Metabolite>("met1", "Metabolite 1");
    auto met2 = std::make_shared<Metabolite>("met2", "Metabolite 2");

    scenario->add_metabolite(met1);
    scenario->add_metabolite(met2);

    const auto& metabolites = scenario->metabolites();
    EXPECT_EQ(metabolites.size(), 2);
    EXPECT_EQ(metabolites[0]->name(), "met1");
    EXPECT_EQ(metabolites[1]->name(), "met2");
}

// Test reactions accessor
TEST_F(ScenarioTest, ReactionsAccessor) {
    auto rxn1 = std::make_shared<Reaction>("R1", 1.0, 1000.0, true, true);
    auto rxn2 = std::make_shared<Reaction>("R2", 2.0, 500.0, true, true);

    scenario->add_reaction(rxn1);
    scenario->add_reaction(rxn2);

    const auto& reactions = scenario->reactions();
    EXPECT_EQ(reactions.size(), 2);
    EXPECT_EQ(reactions[0]->name(), "R1");
    EXPECT_EQ(reactions[1]->name(), "R2");
}

// Test empty scenario
TEST_F(ScenarioTest, EmptyScenario) {
    EXPECT_EQ(scenario->num_metabolites(), 0);
    EXPECT_EQ(scenario->num_reactions(), 0);
    EXPECT_TRUE(scenario->metabolites().empty());
    EXPECT_TRUE(scenario->reactions().empty());
}

// Test adding metabolites and reactions together
TEST_F(ScenarioTest, AddMetabolitesAndReactions) {
    // Add metabolites
    auto glucose = std::make_shared<Metabolite>("glc_D", "D-Glucose");
    auto atp = std::make_shared<Metabolite>("atp", "ATP");
    scenario->add_metabolite(glucose);
    scenario->add_metabolite(atp);

    // Add reactions
    auto rxn = std::make_shared<Reaction>("R1", 1.0, 1000.0, true, true);
    scenario->add_reaction(rxn);

    EXPECT_EQ(scenario->num_metabolites(), 2);
    EXPECT_EQ(scenario->num_reactions(), 1);
}
