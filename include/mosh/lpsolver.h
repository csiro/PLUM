#pragma once

#include "mosh/gapsolver.h"
#include "mosh/lpsolverimp.h"
#include "mosh/dualvals.h"

#include "lime/debug.h"

namespace mosh
{
    class LPSolver : public GapSolver 
    {
    public:
        LPSolver(
            Scenario* scenario, Params* params, int seed,
            int which_formulation, 
            LPSolverImpPtr imp
        ) :
            GapSolver(scenario, params, seed),
            which_experiment_(0), 
            which_formulation_(which_formulation),
            imp_(imp),
            target_flux_(0.0f),
            biomass_search_iters_(0),
            biomass_obj_mult_(params->init_biomass_obj_mult)
        {
        }

        size_t which_experiment() const {return which_experiment_;}
        void set_which_experiment (size_t which_experiment) {
            which_experiment_ = which_experiment;
        }
        size_t which_formulation() const {return which_formulation_;}
        void set_which_formulation (size_t which_formulation);
        
        SolutionPtr solve() override;

        // ALternatively, solve via formulate,
        void formulate (const Experiment* exp);
        // then re-solve after adding reactions
        SolutionPtr do_solve();
        
        // Enable a reaction 
        void enable_reaction (const Reaction* react);
        // Disable a reaction 
        void disable_reaction (const Reaction* react);

        void set_all_react_cost (size_t exp, double cost)
        {
            DEBUG (
                'l', "Set all react cost for exp " << exp << " to " << cost
            );
            for (auto& react_ptr : scenario_->reactions()) {
                Reaction* react = react_ptr.get();
                if (!react->is_biomass())
                    imp_->set_react_cost (exp, react, cost);
            }
        }
        void set_react_cost (size_t exp, Reaction* react, double cost)
        {
            DEBUG (
                'l', "Set react cost for exp " << exp <<
                " react " << react->name() << " to " << cost
            );
            assert (!react->is_biomass());
            imp_->set_react_cost (exp, react, cost);
        }
        
        double biomass_obj_mult() const {return biomass_obj_mult_;}
        // Set obj multiplier. NOTE : Negates the sign
        void set_biomass_obj_mult (size_t exp, double mult)
        {
            DEBUG (
                'l', "Set biomass obj mult for exp " << exp << " to " << mult <<
                " was " << biomass_obj_mult_
            );
            biomass_obj_mult_ = mult;
            imp_->set_biomass_mult (exp, -mult);
        }
        double target_flux() const {return target_flux_;}
        void set_target_flux (double target_flux) {
            target_flux_ = target_flux;
        }

        void write_model (std::string fn)
        {
            imp_->write_model (fn);
        }
        
        std::string summary() override;
        
        void set_quiet (bool quiet) override
        {
            GapSolver::set_quiet (quiet);
            imp_->local_set_quiet (quiet);
        }
        

    protected:
        SolutionPtr biomass_obj_search (SolutionPtr sol);

        size_t which_experiment_;
        size_t which_formulation_;
        LPSolverImpPtr imp_;

        // Target flux for forumation 2
        double target_flux_;
        
        // Solution attributes
        int biomass_search_iters_;
        double biomass_obj_mult_;
        static int num_model_writes_;
    };
    using LPSolverPtr = std::shared_ptr<LPSolver>;
}
