#pragma once

#include <memory>
#include <vector>
#include <iostream>
#include <iomanip>

#include "lime/numutil.h"

#include "mosh/metabolite.h"

namespace mosh
{
    class Experiment
    {
    public:
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
        
        std::string name() const {return name_;}
        size_t index() const {return index_;}
        double biolog_score() const {return biolog_score_;}
        void set_biolog_score (double biolog_score, bool is_growth)
        {
            // Set to 0 if negative
            biolog_score_ = limeMax (0.0, biolog_score);
            is_growth_ = is_growth;
        }
        size_t biolog_rank() const {return biolog_rank_;}
        void set_biolog_rank (size_t biolog_rank) {
            biolog_rank_ = biolog_rank;
        }
        bool is_growth() const {return is_growth_;}
        double base_flux() const {return base_flux_;}
        void set_base_flux (double base_flux) {
            base_flux_ = base_flux;
            has_base_flux_ = true;
        }
        bool has_base_flux() const {return has_base_flux_;}
        
        const std::list<const Metabolite*> carbon_sources() const  {
            return carbon_sources_;
        }
        size_t num_carbon_sources () {
            return carbon_sources_.size();
        }
        void add_carbon_source (const Metabolite* met) {
            carbon_sources_.push_back(met);
        }

        bool is_supply(const Metabolite* met) const
        {
            return is_supply (met->index());
        }
        bool is_supply(size_t met_idx) const
        {
            return limeIsNegative(lb_[met_idx]);
        }
        double supply(const Metabolite* met) const
        {
            return supply (met->index());
        }
        double supply(size_t met_idx) const
        {
            assert (is_supply(met_idx));
            return -lb_[met_idx];
        }
        
        bool is_residual(const Metabolite* met) const
        {
            return is_residual(met->index());
        }
        bool is_residual(size_t met_idx) const
        {
            return limeIsPositive(ub_[met_idx]);
        }
        double residual(const Metabolite* met) const
        {
            return residual(met->index());
        }
        double residual(size_t met_idx) const
        {
            assert (is_residual(met_idx));
            return ub_[met_idx];
        }
        bool is_balanced(const Metabolite* met) const
        {
            return is_balanced(met->index());
        }
        bool is_balanced(size_t met_idx) const
        {
            return
                limeIsZero(lb_[met_idx]) &&
                limeIsZero(ub_[met_idx]);
        }

        double lb (const Metabolite* met) const 
        {
            return lb (met->index());
        }
        double lb (size_t met_idx) const
        {
            return lb_[met_idx];
        }
        double ub (const Metabolite* met) const
        {
            return ub (met->index());
        }
        double ub (size_t met_idx) const
        {
            return ub_[met_idx];
        }
        
        // When set, this value has human-centred meaning:
        // +ve -> this is a supplied metabolite
        // -ve -> this is a residual metabolite
        //   0 -> this is a balanced metabolite
        void set_supply_resid (const Metabolite* met, double lb, double ub)
        {
            set_supply_resid (met->index(), lb, ub);
        }
        void set_supply_resid (size_t met_idx, double lb, double ub)
        {
            // Negate to save in the mathematically-centred meaning
            lb_[met_idx] = limeMin (-lb, -ub);
            ub_[met_idx] = limeMax (-lb, -ub);
        }

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
        
        std::string name_;
        size_t index_;
        double biolog_score_;
        // In the experiments, where do I sit?
        // biolog_rank_ == 0 -> I have the maximum biolog score
        size_t biolog_rank_;
        double base_flux_;
        bool has_base_flux_;
        bool is_growth_;
        std::list<const Metabolite*> carbon_sources_;
        std::vector<double> lb_;
        std::vector<double> ub_;
    };
    
    using ExperimentPtr = std::shared_ptr<Experiment>;
}
