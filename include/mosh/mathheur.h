#pragma once

// Compile with
// g++ fn.cpp -lgurobi81 -lgurobi_g++5.2

#include <vector>
#include <memory>

#include "gurobi_c++.h"

#include "mosh/gapsolver.h"
#include "mosh/grblpsolver.h"

#include "lime/numutil.h"
#include "lime/rand.h"

namespace mosh
{
    class MathHeur : public GapSolver 
    {
    public:
        MathHeur (
            Scenario* scenario, Params* params, 
            GrbEnvPtr env, int seed
        ) :
            GapSolver(scenario, params, seed),
            iter_(0),
            max_iters_(params->max_mh_iters),
            model_(*env),
            env_(env),
            use_(num_reactions()),
            flux_(num_reactions()),
            biomass_()
        {
            model_.set(GRB_StringAttr_ModelName, "plum_mh");
            model_.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
            model_.set(GRB_IntParam_Threads, params->num_threads);
            model_.set(GRB_IntParam_Seed, grb_seed());
            model_.set(GRB_IntParam_MIPFocus, GRB_MIPFOCUS_FEASIBILITY);
            model_.set(GRB_DoubleParam_MIPGap, params->max_mip_gap);
            model_.set(GRB_DoubleParam_FeasibilityTol, params->feasibility_tol);
            model_.set(GRB_DoubleParam_IntFeasTol, params->int_feas_tol);
            if (params->time_limit > 0) {
                double dbl_limit = (double) params->time_limit;
                model_.set(GRB_DoubleParam_TimeLimit, dbl_limit);
            }
            //limeSetEpsilon(1e-8);
        }

        SolutionPtr solve() override;
        
        GRBVar& flux(size_t k) {return flux_[k];}
        
        std::string summary() override;
        
        void write_model (std::string filename)
        {
            model_.write (filename);
        }

    private:
        void set_up_model();
        void improve(const Metabolite* met, SolutionPtr& sol);
        
        int iter_;
        int max_iters_;
        int num_threads_;
        int time_limit_;
        lime::Rand rand_;
        
        GRBModel model_;
        GrbEnvPtr env_;
        
        // Binary var for Do I use this rection?
        std::vector<GRBVar> use_;
        // Variable for flux through a reaction
        std::vector<GRBVar> flux_;

        // The biomass var
        GRBVar biomass_;

    };
}
