/**
 * @file reaction_test.cpp
 * @brief Unit tests for the Reaction class using Google Test
 */

#include <gtest/gtest.h>
#include <sstream>
#include "mosh/reaction.h"
#include "mosh/metabolite.h"

using namespace mosh;

/**
 * @class ReactionTest
 * @brief Test fixture for Reaction class tests
 */
class ReactionTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create test metabolites
        glucose = new Metabolite("glc_D", "D-Glucose");
        glucose->set_index(0);

        atp = new Metabolite("atp", "ATP");
        atp->set_index(1);

        adp = new Metabolite("adp", "ADP");
        adp->set_index(2);

        h2o = new Metabolite("h2o", "Water");
        h2o->set_index(3);

        // Create test reactions
        biomass = new Reaction("BIOMASS", -1.0, 1000.0, true, true, false);
        exchange = new Reaction("EX_glc", 0.0, 10.0, true, true, false);
        dummy = new Reaction("DUMMY_R", 0.0, 1000.0, false, false, true);
    }

    void TearDown() override {
        delete glucose;
        delete atp;
        delete adp;
        delete h2o;
        delete biomass;
        delete exchange;
        delete dummy;
    }

    Metabolite* glucose;
    Metabolite* atp;
    Metabolite* adp;
    Metabolite* h2o;
    Reaction* biomass;
    Reaction* exchange;
    Reaction* dummy;
};

// Test constructor and basic accessors
TEST_F(ReactionTest, ConstructorInitializesCorrectly) {
    EXPECT_EQ(biomass->name(), "BIOMASS");
    EXPECT_EQ(biomass->obj_coeff(), -1.0);
    EXPECT_EQ(biomass->flux_ub(), 1000.0);
    EXPECT_TRUE(biomass->is_active());
    EXPECT_TRUE(biomass->is_selected());
    EXPECT_FALSE(biomass->is_dummy());
    EXPECT_TRUE(biomass->is_biomass());  // negative obj_coeff
}

// Test name setters and getters
TEST_F(ReactionTest, NameSetterGetter) {
    biomass->set_name("BIOMASS_2");
    EXPECT_EQ(biomass->name(), "BIOMASS_2");

    biomass->set_full_name("Biomass Reaction Version 2");
    EXPECT_EQ(biomass->full_name(), "Biomass Reaction Version 2");
}

// Test index setter and getter
TEST_F(ReactionTest, IndexSetterGetter) {
    biomass->set_index(10);
    EXPECT_EQ(biomass->index(), 10);
    EXPECT_EQ(biomass->lps_id(), 11);  // 1-based for LP solver
}

// Test objective coefficient
TEST_F(ReactionTest, ObjectiveCoefficientSetterGetter) {
    EXPECT_EQ(biomass->obj_coeff(), -1.0);
    biomass->set_obj_coeff(-2.5);
    EXPECT_EQ(biomass->obj_coeff(), -2.5);
}

// Test biomass flag based on obj_coeff
TEST_F(ReactionTest, BiomassFlag) {
    EXPECT_TRUE(biomass->is_biomass());  // obj_coeff = -1.0
    EXPECT_FALSE(exchange->is_biomass()); // obj_coeff = 0.0
}

// Test active flag
TEST_F(ReactionTest, ActiveFlagSetterGetter) {
    EXPECT_TRUE(biomass->is_active());
    biomass->set_active(false);
    EXPECT_FALSE(biomass->is_active());
}

// Test selected flag
TEST_F(ReactionTest, SelectedFlagSetterGetter) {
    EXPECT_TRUE(biomass->is_selected());
    biomass->set_selected(false);
    EXPECT_FALSE(biomass->is_selected());
}

// Test dummy flag
TEST_F(ReactionTest, DummyFlag) {
    EXPECT_FALSE(biomass->is_dummy());
    EXPECT_TRUE(dummy->is_dummy());
}

// Test flux upper bound
TEST_F(ReactionTest, FluxUpperBound) {
    EXPECT_EQ(biomass->flux_ub(), 1000.0);
    biomass->set_flux_ub(500.0);
    EXPECT_EQ(biomass->flux_ub(), 500.0);
}

// Test known flux
TEST_F(ReactionTest, KnownFlux) {
    EXPECT_FALSE(biomass->has_known_flux());

    biomass->set_known_flux(5.5, 0.1);
    EXPECT_TRUE(biomass->has_known_flux());
    EXPECT_EQ(biomass->known_flux(), 5.5);
    EXPECT_EQ(biomass->known_flux_error(), 0.1);
}

// Test reduced cost
TEST_F(ReactionTest, ReducedCost) {
    EXPECT_EQ(biomass->reduced_cost(), 0.0);
    biomass->set_reduced_cost(3.14);
    EXPECT_EQ(biomass->reduced_cost(), 3.14);
}

// Test setting metabolite coefficients
TEST_F(ReactionTest, SetCoefficient) {
    // ATP + H2O -> ADP
    exchange->set_coeff(atp, -1.0);
    exchange->set_coeff(h2o, -1.0);
    exchange->set_coeff(adp, 1.0);

    EXPECT_EQ(exchange->num_mets(), 3);
    EXPECT_TRUE(exchange->uses(atp));
    EXPECT_TRUE(exchange->uses(h2o));
    EXPECT_TRUE(exchange->uses(adp));
    EXPECT_FALSE(exchange->uses(glucose));

    EXPECT_EQ(exchange->met_coeff(atp), -1.0);
    EXPECT_EQ(exchange->met_coeff(h2o), -1.0);
    EXPECT_EQ(exchange->met_coeff(adp), 1.0);
}

// Test in_mets and out_mets
TEST_F(ReactionTest, InputOutputMetabolites) {
    // ATP + H2O -> ADP
    exchange->set_coeff(atp, -1.0);
    exchange->set_coeff(h2o, -1.0);
    exchange->set_coeff(adp, 1.0);

    EXPECT_EQ(exchange->in_mets().size(), 2);
    EXPECT_EQ(exchange->out_mets().size(), 1);
}

// Test EX reaction (no reactants)
TEST_F(ReactionTest, ExchangeReaction) {
    Reaction ex_reaction("EX_glc", 0.0, 10.0, true, true);
    ex_reaction.set_coeff(glucose, 1.0);  // Only products

    EXPECT_TRUE(ex_reaction.is_ex());
    EXPECT_FALSE(ex_reaction.is_dm());
    EXPECT_TRUE(ex_reaction.is_ex_dm());
}

// Test DM reaction (no products)
TEST_F(ReactionTest, DemandReaction) {
    Reaction dm_reaction("DM_atp", 0.0, 10.0, true, true);
    dm_reaction.set_coeff(atp, -1.0);  // Only reactants

    EXPECT_FALSE(dm_reaction.is_ex());
    EXPECT_TRUE(dm_reaction.is_dm());
    EXPECT_TRUE(dm_reaction.is_ex_dm());
}

// Test export reaction flag (all positive coefficients)
TEST_F(ReactionTest, ExportReaction) {
    Reaction export_rxn("EXPORT", 0.0, 10.0, true, true);
    export_rxn.set_coeff(glucose, 1.0);
    export_rxn.set_coeff(atp, 2.0);

    EXPECT_TRUE(export_rxn.is_export());
}

// Test non-export reaction (has negative coefficients)
TEST_F(ReactionTest, NonExportReaction) {
    Reaction normal_rxn("NORMAL", 0.0, 10.0, true, true);
    normal_rxn.set_coeff(atp, -1.0);
    normal_rxn.set_coeff(adp, 1.0);

    EXPECT_FALSE(normal_rxn.is_export());
}

// Test formula generation
TEST_F(ReactionTest, FormulaGeneration) {
    Reaction rxn("R1", 0.0, 10.0, true, true);
    rxn.set_coeff(atp, -1.0);
    rxn.set_coeff(h2o, -1.0);
    rxn.set_coeff(adp, 1.0);

    std::string formula = rxn.formula();
    EXPECT_TRUE(formula.find("atp") != std::string::npos);
    EXPECT_TRUE(formula.find("h2o") != std::string::npos);
    EXPECT_TRUE(formula.find("adp") != std::string::npos);
    EXPECT_TRUE(formula.find("->") != std::string::npos);
}

// Test met_has_coeff
TEST_F(ReactionTest, MetHasCoeff) {
    exchange->set_coeff(glucose, -1.0);

    EXPECT_TRUE(exchange->met_has_coeff(glucose->index()));
    EXPECT_FALSE(exchange->met_has_coeff(atp->index()));
}

// Test coefficient by index
TEST_F(ReactionTest, CoefficientByIndex) {
    exchange->set_coeff(glucose, -2.5);

    EXPECT_EQ(exchange->met_coeff(glucose->index()), -2.5);
    EXPECT_EQ(exchange->met_coeff(999), 0.0);  // Non-existent metabolite
}

// Test remove_met
TEST_F(ReactionTest, RemoveMetabolite) {
    exchange->set_coeff(atp, -1.0);
    exchange->set_coeff(adp, 1.0);
    EXPECT_EQ(exchange->num_mets(), 2);

    exchange->remove_met(atp);
    EXPECT_EQ(exchange->num_mets(), 2);  // Still has entry, but coeff = 0
    EXPECT_FALSE(exchange->uses(atp));
    EXPECT_EQ(exchange->met_coeff(atp), 0.0);
}

// Test same_as comparison
TEST_F(ReactionTest, SameAsComparison) {
    Reaction rxn1("R1", 0.0, 10.0, true, true);
    rxn1.set_coeff(atp, -1.0);
    rxn1.set_coeff(adp, 1.0);

    Reaction rxn2("R2", 0.0, 10.0, true, true);
    rxn2.set_coeff(atp, -1.0);
    rxn2.set_coeff(adp, 1.0);

    EXPECT_TRUE(rxn1.same_as(&rxn2));
}

// Test same_as with different reactions
TEST_F(ReactionTest, SameAsDifferent) {
    Reaction rxn1("R1", 0.0, 10.0, true, true);
    rxn1.set_coeff(atp, -1.0);
    rxn1.set_coeff(adp, 1.0);

    Reaction rxn2("R2", 0.0, 10.0, true, true);
    rxn2.set_coeff(atp, -2.0);  // Different coefficient
    rxn2.set_coeff(adp, 1.0);

    EXPECT_FALSE(rxn1.same_as(&rxn2));
}

// Test reverse_of comparison
TEST_F(ReactionTest, ReverseOfComparison) {
    Reaction rxn1("R1", 0.0, 10.0, true, true);
    rxn1.set_coeff(atp, -1.0);
    rxn1.set_coeff(adp, 1.0);

    Reaction rxn2("R2", 0.0, 10.0, true, true);
    rxn2.set_coeff(atp, 1.0);   // Opposite sign
    rxn2.set_coeff(adp, -1.0);  // Opposite sign

    EXPECT_TRUE(rxn1.reverse_of(&rxn2));
    EXPECT_TRUE(rxn2.reverse_of(&rxn1));
}

// Test write_to output
TEST_F(ReactionTest, WriteToOutputStream) {
    Reaction rxn("R1", 5.0, 100.0, true, true);
    rxn.set_full_name("Test Reaction 1");
    rxn.set_coeff(atp, -1.0);
    rxn.set_coeff(adp, 1.0);

    std::ostringstream oss;
    rxn.write_to(oss);

    std::string output = oss.str();
    EXPECT_TRUE(output.find("REACTION R1") != std::string::npos);
    EXPECT_TRUE(output.find("5") != std::string::npos);  // obj_coeff
    EXPECT_TRUE(output.find("100") != std::string::npos);  // flux_ub
}

// Test nicify function
TEST_F(ReactionTest, NicifyFunction) {
    Reaction rxn("R1", 0.0, 10.0, true, true);

    Metabolite met_encoded("compound__40__test__41__", "Test");
    met_encoded.set_index(10);
    rxn.set_coeff(&met_encoded, 1.0);

    std::string nice = rxn.nice_formula();
    EXPECT_TRUE(nice.find("(") != std::string::npos);
    EXPECT_TRUE(nice.find(")") != std::string::npos);
}
