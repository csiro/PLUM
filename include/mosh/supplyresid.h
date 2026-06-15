/**
 * @file supplyresid.h
 * @brief Supply and residual metabolite representation for FBA gap-filling
 *
 * This file defines the SupplyResid class which tracks metabolite supply/residual
 * bounds in the context of metabolic network gap-filling. Metabolites can be
 * supplied (available as inputs), residual (produced as outputs), or balanced
 * (internally conserved) in the flux balance analysis.
 */
#pragma once

#include <memory>

/**
 * @namespace mosh
 * @brief Main namespace for PLUM metabolic optimization solver components
 */
namespace mosh
{
/**
 * @class SupplyResid
 * @brief Represents supply and residual bounds for a metabolite in FBA
 *
 * Tracks whether a metabolite is supplied (available as input), residual
 * (produced as output), or balanced (mass-balanced internally). Uses internal
 * mathematical representation where negative values indicate supply and positive
 * values indicate residual, with conversions to human-readable semantics.
 *
 * The class maintains both lower and upper bounds for supply/residual constraints
 * in the linear programming formulation of the metabolic network.
 */
    class SupplyResid 
    {
    public:
    /**
     * @brief Constructs a SupplyResid object with initial bounds
     * @param lb Lower bound for the metabolite constraint
     * @param ub Upper bound for the metabolite constraint
     */
        SupplyResid(double lb, double ub) :
            lb_(lb),
            ub_(ub)
        {
        }

        /**
         * @brief Gets the index of this metabolite in the problem formulation
         * @return The zero-based index of this metabolite
         */
        size_t index() const {return index_;}
        /**
         * @brief Sets the index of this metabolite in the problem formulation
         * @param index The zero-based index to assign
         */
        void set_index (size_t index) {index_ = index;}
        /**
         * @brief Checks if this metabolite is supplied (available as input)
         * @return true if the metabolite has a supply constraint (negative lower bound)
         */
        bool is_supply() const {return limeIsNegative(sr_lb_);}
        /**
         * @brief Gets the supply amount for this metabolite
         * @return The supply quantity (positive value)
         * @pre is_supply() must return true
         */
        double supply() const {
            assert (is_supply());
            return -sr_lb_;
        }
        /**
         * @brief Checks if this metabolite is residual (produced as output)
         * @return true if the metabolite has a residual constraint (positive upper bound)
         */
        bool is_residual() const {return limeIsPositive(sr_ub_);}
        /**
         * @brief Gets the residual amount for this metabolite
         * @return The residual quantity (positive value)
         * @pre is_residual() must return true
         */
        double residual() const {
            assert (is_residual());
            return sr_ub_;
        }
        /**
         * @brief Checks if this metabolite is mass-balanced (neither supplied nor residual)
         * @return true if the metabolite has zero supply/residual bounds
         */
        bool is_balanced() const {return limeIsZero(sr_lb_);}
        /**
         * @brief Gets the lower bound in mathematical representation
         * @return The lower bound (negative indicates supply)
         */
        double sr_lb() const {return sr_lb_;}
        /**
         * @brief Gets the upper bound in mathematical representation
         * @return The upper bound (positive indicates residual)
         */
        double sr_ub() const {return sr_ub_;}
        /**
         * @brief Sets supply/residual bounds using human-readable convention
         * @param lb Lower bound in human convention (+ve = supply, -ve = residual, 0 = balanced)
         * @param ub Upper bound in human convention (+ve = supply, -ve = residual, 0 = balanced)
         *
         * Converts from human-centered convention (positive = supply) to mathematical
         * convention (negative = supply) used internally for LP formulation.
         */
        void set_supply_resid (double lb, double ub) {
            // Negate to save in the mathematically-centred meaning
            sr_lb_ = limeMin (-lb, -ub);
            sr_ub_ = limeMax (-lb, -ub);
            DEBUG (
                'm', "  Set lb/ub to " << sr_lb_ << " " << sr_ub_ <<
                " " << sr_ub() 
            );
        }
    
    private:
        std::string name_; /**< Name identifier for this metabolite */"
        size_t index_; /**< Index of this metabolite in the LP problem formulation */"
        /**
         * @brief Supply/residual lower bound (mathematical representation)
         *
         * Internal convention:
         * - Negative value: supplied metabolite
         * - Positive value: residual metabolite
         * - Zero: balanced metabolite
         */
        double sr_lb_;  
        /**
         * @brief Supply/residual upper bound (mathematical representation)
         *
         * Internal convention:
         * - Negative value: supplied metabolite
         * - Positive value: residual metabolite
         * - Zero: balanced metabolite
         */
        double sr_ub_;
        double dummy_; /**< Dummy variable for LP formulation purposes */"
    };
    
    /**
     * @typedef SupplyResidPtr
     * @brief Shared pointer type for SupplyResid objects
     *
     * Used for managing lifetime and sharing of SupplyResid instances
     * across the metabolic network representation.
     */
    using SupplyResidPtr = std::shared_ptr<SupplyResid>;
}
