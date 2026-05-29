#pragma once

/*
  An incremental solver.
  Adds reactions into the model in order of increasing cost (obj-coeff)
  until biomass is produced (without dummies)
 */

#include <set>

#include "mosh/gapsolver.h"
#include "mosh/lpsolver.h"
#include "mosh/reaction.h"

namespace mosh
{
    class IncrSolver : public GapSolver 
    {
    public:
        /** Create incremental solver.
            imp2 is non-null only if integer solution required
        */
        IncrSolver(
            Scenario* scenario, Params* params, int seed,
            LPSolverImpPtr imp, LPSolverImpPtr imp2,
            double min_cost, std::ofstream& selected_react_file
        ) :
            GapSolver(scenario, params, seed),
            lp_solver_ (scenario, params, rand_seed(), 1, imp),
            imp2_(imp2),
            return_int_sol_(imp2 != nullptr),
            time_limit_(params->time_limit),
            min_cost_(min_cost),
            max_used_cost_(min_cost),
            costs_(),
            selected_react_file_(selected_react_file),
            available_(scenario->num_metabolites(), false),
            summary_(),
            cost_iters_(0),
            level_iters_(0),
            cycle_mets_(0)
        {
            init();
        }

        SolutionPtr solve() override;
        
        std::string summary() override;
        
    private:
        void init();
        SolutionPtr do_solve();
        // Process given cost
        SolutionPtr cost_iter (double cost);
        // Next level. Return false if no reactions added
        bool level_iter (double cost, int level);
        bool add_cycle_met (double cost);

        bool is_enabled (const Reaction* react);
        bool makes_carbon (const Reaction* react);
        void propagate (const Reaction* react);
        bool is_available (const Metabolite* met) const
        {
            return available_[met->index()];
        }
            

        LPSolver lp_solver_;
        LPSolverImpPtr imp2_;
        bool return_int_sol_;
        int time_limit_;
        double min_cost_;
        double max_used_cost_;
        std::set<double> costs_;
        lime::Rand rand_;
        std::ofstream& selected_react_file_;

        // Array over metabolites
        // Has the metabolite been produced by a reaction at
        // a previous level
        std::vector<bool> available_;

        // Summary stats
        std::string summary_;
        size_t cost_iters_;
        size_t level_iters_;
        size_t cycle_mets_;
    };
}
