/**
 * @file experiment.h
 * @brief Defines the Experiment class for representing metabolic growth experiments
 *
 * This file contains the Experiment class which models individual biological experiments
 * with metabolite supply/residual constraints for flux balance analysis and gap-filling.
 */

#pragma once

#include <memory>
#include <vector>
#include <iostream>
#include <iomanip>

#include "lime/numutil.h"

#include "mosh/metabolite.h"

/**
 * @brief Namespace for metabolic optimization and gap-filling components
 */
namespace mosh
{
/**
 * @class Experiment
 * @brief Represents a single biological growth experiment with metabolite constraints
 *
 * This class encapsulates experimental data including metabolite supply/residual bounds,
 * carbon sources, Biolog phenotype microarray scores, and base flux values. It supports
 * flux balance analysis by tracking metabolite availability constraints (lower/upper bounds)
 * that distinguish between supplied metabolites (negative flux), residual metabolites
 * (positive flux), and balanced metabolites (zero flux).
 */
    class Experiment
    {
    public:
    /**
     * @brief Constructs an Experiment with the given name, index, and metabolite count
     * @param name Name of the experiment
     * @param index Unique index identifying this experiment
     * @param num_metabolites Number of metabolites in the metabolic network
     */
        Experiment(std::string name, size_t index, size_t num_metabolites) :
            name_(name),
            index_(index),
            biolog_score_(1.0f),
            biolog_rank_(0),
            base_flux_(0.0f),
            has_base_flux_(false),
            is_growth_(true),
            carbon_sources_(),
            lb_(num_metabolites),
            ub_(num_metabolites)
        {
        }
        
    /**
     * @brief Gets the experiment name
     * @return Name of the experiment
     */
        std::string name() const {return name_;}
    /**
     * @brief Gets the experiment index
     * @return Unique index of the experiment
     */
        size_t index() const {return index_;}
    /**
     * @brief Gets the Biolog phenotype score
     * @return Biolog score indicating growth intensity (non-negative)
     */
        double biolog_score() const {return biolog_score_;}
    /**
     * @brief Sets the Biolog phenotype score and growth status
     * @param biolog_score Raw Biolog score (clamped to non-negative values)
     * @param is_growth True if the experiment indicates growth, false otherwise
     */
        void set_biolog_score (double biolog_score, bool is_growth)
        {
            // Set to 0 if negative
            biolog_score_ = limeMax (0.0, biolog_score);
            is_growth_ = is_growth;
        }
    /**
     * @brief Gets the Biolog rank among all experiments
     * @return Rank where 0 indicates maximum Biolog score
     */
        size_t biolog_rank() const {return biolog_rank_;}
    /**
     * @brief Sets the Biolog rank
     * @param biolog_rank Rank of this experiment (0 = highest score)
     */
        void set_biolog_rank (size_t biolog_rank) {
            biolog_rank_ = biolog_rank;
        }
    /**
     * @brief Checks if this experiment indicates growth
     * @return True if growth was observed, false otherwise
     */
        bool is_growth() const {return is_growth_;}
    /**
     * @brief Gets the base flux value
     * @return Base flux (e.g., biomass flux) for this experiment
     */
        double base_flux() const {return base_flux_;}
    /**
     * @brief Sets the base flux value
     * @param base_flux Base flux value (e.g., predicted biomass production rate)
     */
        void set_base_flux (double base_flux) {
            base_flux_ = base_flux;
            has_base_flux_ = true;
        }
    /**
     * @brief Checks if a base flux has been set
     * @return True if base flux has been assigned, false otherwise
     */
        bool has_base_flux() const {return has_base_flux_;}
        
    /**
     * @brief Gets the list of carbon sources for this experiment
     * @return List of metabolite pointers representing carbon sources
     */
        const std::list<const Metabolite*> carbon_sources() const  {
            return carbon_sources_;
        }
    /**
     * @brief Gets the number of carbon sources
     * @return Count of carbon source metabolites
     */
        size_t num_carbon_sources () {
            return carbon_sources_.size();
        }
    /**
     * @brief Adds a carbon source metabolite to this experiment
     * @param met Pointer to the carbon source metabolite
     */
        void add_carbon_source (const Metabolite* met) {
            carbon_sources_.push_back(met);
        }

    /**
     * @brief Checks if a metabolite is supplied in this experiment
     * @param met Pointer to the metabolite to check
     * @return True if the metabolite is supplied (negative lower bound)
     */
        bool is_supply(const Metabolite* met) const
        {
            return is_supply (met->index());
        }
    /**
     * @brief Checks if a metabolite is supplied in this experiment
     * @param met_idx Index of the metabolite to check
     * @return True if the metabolite is supplied (negative lower bound)
     */
        bool is_supply(size_t met_idx) const
        {
            return limeIsNegative(lb_[met_idx]);
        }
    /**
     * @brief Gets the supply amount for a metabolite
     * @param met Pointer to the supplied metabolite
     * @return Supply amount (positive value)
     */
        double supply(const Metabolite* met) const
        {
            return supply (met->index());
        }
    /**
     * @brief Gets the supply amount for a metabolite
     * @param met_idx Index of the supplied metabolite
     * @return Supply amount (positive value)
     */
        double supply(size_t met_idx) const
        {
            assert (is_supply(met_idx));
            return -lb_[met_idx];
        }
        
    /**
     * @brief Checks if a metabolite is residual (produced/secreted) in this experiment
     * @param met Pointer to the metabolite to check
     * @return True if the metabolite is residual (positive upper bound)
     */
        bool is_residual(const Metabolite* met) const
        {
            return is_residual(met->index());
        }
    /**
     * @brief Checks if a metabolite is residual (produced/secreted) in this experiment
     * @param met_idx Index of the metabolite to check
     * @return True if the metabolite is residual (positive upper bound)
     */
        bool is_residual(size_t met_idx) const
        {
            return limeIsPositive(ub_[met_idx]);
        }
    /**
     * @brief Gets the residual amount for a metabolite
     * @param met Pointer to the residual metabolite
     * @return Residual amount (positive value)
     */
        double residual(const Metabolite* met) const
        {
            return residual(met->index());
        }
    /**
     * @brief Gets the residual amount for a metabolite
     * @param met_idx Index of the residual metabolite
     * @return Residual amount (positive value)
     */
        double residual(size_t met_idx) const
        {
            assert (is_residual(met_idx));
            return ub_[met_idx];
        }
    /**
     * @brief Checks if a metabolite is balanced (neither supplied nor residual)
     * @param met Pointer to the metabolite to check
     * @return True if the metabolite has zero bounds (balanced)
     */
        bool is_balanced(const Metabolite* met) const
        {
            return is_balanced(met->index());
        }
    /**
     * @brief Checks if a metabolite is balanced (neither supplied nor residual)
     * @param met_idx Index of the metabolite to check
     * @return True if the metabolite has zero bounds (balanced)
     */
        bool is_balanced(size_t met_idx) const
        {
            return
                limeIsZero(lb_[met_idx]) &&
                limeIsZero(ub_[met_idx]);
        }

    /**
     * @brief Gets the lower bound for a metabolite
     * @param met Pointer to the metabolite
     * @return Lower bound (negative indicates supply)
     */
        double lb (const Metabolite* met) const 
        {
            return lb (met->index());
        }
    /**
     * @brief Gets the lower bound for a metabolite
     * @param met_idx Index of the metabolite
     * @return Lower bound (negative indicates supply)
     */
        double lb (size_t met_idx) const
        {
            return lb_[met_idx];
        }
    /**
     * @brief Gets the upper bound for a metabolite
     * @param met Pointer to the metabolite
     * @return Upper bound (positive indicates residual)
     */
        double ub (const Metabolite* met) const
        {
            return ub (met->index());
        }
    /**
     * @brief Gets the upper bound for a metabolite
     * @param met_idx Index of the metabolite
     * @return Upper bound (positive indicates residual)
     */
        double ub (size_t met_idx) const
        {
            return ub_[met_idx];
        }
        
        // When set, this value has human-centred meaning:
        // +ve -> this is a supplied metabolite
        // -ve -> this is a residual metabolite
        //   0 -> this is a balanced metabolite
    /**
     * @brief Sets supply/residual bounds for a metabolite using human-centered values
     * @param met Pointer to the metabolite
     * @param lb Lower bound (positive = supply, negative = residual, 0 = balanced)
     * @param ub Upper bound (positive = supply, negative = residual, 0 = balanced)
     *
     * Note: The method converts human-centered values (positive = supply) to
     * mathematical representation (negative = supply) for internal storage.
     */
        void set_supply_resid (const Metabolite* met, double lb, double ub)
        {
            set_supply_resid (met->index(), lb, ub);
        }
    /**
     * @brief Sets supply/residual bounds for a metabolite using human-centered values
     * @param met_idx Index of the metabolite
     * @param lb Lower bound (positive = supply, negative = residual, 0 = balanced)
     * @param ub Upper bound (positive = supply, negative = residual, 0 = balanced)
     *
     * Note: The method converts human-centered values (positive = supply) to
     * mathematical representation (negative = supply) for internal storage.
     */
        void set_supply_resid (size_t met_idx, double lb, double ub)
        {
            // Negate to save in the mathematically-centred meaning
            lb_[met_idx] = limeMin (-lb, -ub);
            ub_[met_idx] = limeMax (-lb, -ub);
        }

    /**
     * @brief Displays experiment details to an output stream
     * @param out Output stream to write to
     */
        void show (std::ostream& out)
        {
            out << "Experiment " << name_ << " (" << index_ << ")" << std::endl;
            for (size_t k = 0; k < lb_.size(); k++) {
                if (!is_balanced(k)) {
                    out << "  Met " << std::setw(3) << k <<
                        " [" << lb(k) << "," << ub(k) << "]" << std::endl;
                }
            }
        }
    
    private:
        // Store lower and upper bounds
        // -ve -> this is a supplied metabolite
        // +ve -> this is a residual metabolite
        //   0 -> this is a balanced metabolite
        
        std::string name_; /**< Name of the experiment **/
        size_t index_; /**< Unique index identifying this experiment **/
        double biolog_score_; /**< Biolog phenotype microarray score (non-negative) **/
        // In the experiments, where do I sit?
        // biolog_rank_ == 0 -> I have the maximum biolog score
        size_t biolog_rank_; /**< Rank among experiments (0 = highest Biolog score) **/
        double base_flux_; /**< Base flux value (e.g., biomass production rate) **/
        bool has_base_flux_; /**< Flag indicating if base flux has been set **/
        bool is_growth_; /**< True if the experiment indicates growth **/
        std::list<const Metabolite*> carbon_sources_; /**< List of carbon source metabolites **/
        std::vector<double> lb_; /**< Lower bounds for metabolites (negative = supply) **/
        std::vector<double> ub_; /**< Upper bounds for metabolites (positive = residual) **/
    };
    
    /**
     * @typedef ExperimentPtr
     * @brief Shared pointer type for Experiment objects
     */
    using ExperimentPtr = std::shared_ptr<Experiment>;
}
