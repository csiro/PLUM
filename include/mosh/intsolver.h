#pragma once

#include <vector>
#include <memory>

#include "mosh/gapsolver.h"
#include "mosh/lpsolverimp.h"

namespace mosh
{
    class IntSolver : public GapSolver 
    {
    public:
        IntSolver(
            Scenario* scenario, Params* params, int seed,
            int which_formulation, 
            LPSolverImpPtr imp, MultiSolPtr cts_sol
        ) :
            GapSolver(scenario, params, seed),
            which_formulation_(which_formulation),
            imp_(imp),
            cts_sol_(cts_sol),
            biomass_obj_mult_(0.0f),
            rounding_iters_(0),
            mip_gap_(-1.0f)
        {
        }

        SolutionPtr solve() override;
        
        std::string summary() override;
        
    private:
        /** Given a solution, find the biomass mult to use for integer
            formulation
         */
        void calc_biomass_mult();
        void formulate ();
        void fix_rounding();

        /* Formulation 1 uses biomass mult; formulation 2 doesn't */
        int which_formulation_;
        LPSolverImpPtr imp_;
        
        /* cts_sol is the solution to to the CTS prob for each experiment
        */
        MultiSolPtr cts_sol_;
        
        // Solution attributes
        double biomass_obj_mult_;
        int rounding_iters_;
        double mip_gap_;
    };
}
