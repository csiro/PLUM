#pragma once

// Compile with
// g++ fn.cpp -lgurobi81 -lgurobi_g++5.2

#include <vector>
#include <memory>

#include "gurobi_c++.h"

#include "mosh/lpsolverimp.h"
#include "mosh/dualvals.h"

#include "lime/numutil.h"

namespace mosh
{
    using GrbEnvPtr = std::shared_ptr<GRBEnv>;
    
    class GrbLPSolverImp : public LPSolverImp
    {
    public:
        GrbLPSolverImp(
            Scenario* scenario, Params* params, int seed, GrbEnvPtr env
        ) :
            LPSolverImp(scenario, params),
            model_(*env),
            cts_sol_(nullptr),
            flux_(
                scenario_->num_reactions(),
                std::vector<GRBVar>(scenario->num_experiments())
            ),
            use_(scenario_->num_reactions()),
            constr_(),
            biomass_(scenario->num_experiments())
        {
            model_.set(GRB_StringAttr_ModelName, "plum_lp");
            model_.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
            model_.set(GRB_IntParam_Threads, params->num_threads);
            model_.set(GRB_IntParam_Seed, seed);
            if (params->time_limit > 0) {
                double dbl_limit = (double) params->time_limit;
                model_.set(GRB_DoubleParam_TimeLimit, dbl_limit);
            }
            model_.set(GRB_DoubleParam_MIPGap, params_->max_mip_gap);
            model_.set(GRB_DoubleParam_IntFeasTol, params_->int_feas_tol);
            model_.set(GRB_DoubleParam_FeasibilityTol, params_->feasibility_tol);
            //limeSetEpsilon(1e-8);
        }

        void init_cts() override
        {
        }
        
        void init_int(SolutionPtr cts_sol) override
        {
            cts_sol_ = cts_sol;
            // Switch off presolve... (for debugging)
            //model_.set(GRB_IntParam_Presolve, 0);
        }

        void make_flux_var (
            const Reaction* react, size_t exp, std::string name,
            double obj_coeff, double lb, double ub
        ) override;
        void make_use_var (
            const Reaction* react, std::string name, double obj_coeff, double lb, double ub
        ) override;
        void add_met_constraint (
            const Metabolite* met, size_t experiment, double lb, double ub,
            std::vector<Reaction*>& reacts, std::vector<double>& coeffs
        ) override;
        void add_react_link_constraint (
            const Reaction* react, size_t exp
        ) override;
        void set_react_bounds (
            const Reaction* react, size_t exp, double lb, double ub
        ) override
        {
            DEBUG (
                'z', "      Set bounds for react " << *react <<
                " exp " << exp <<
                " to [" << lb << "," << ub << "]"
            );
            flux_[react->index()][exp].set (GRB_DoubleAttr_LB, lb);
            flux_[react->index()][exp].set (GRB_DoubleAttr_UB, ub);
        }
        void set_react_cost (
            size_t exp, const Reaction* react, double cost
        ) override
        {
            DEBUG (
                'z', "      Set cost for react " << *react <<
                " to " << cost
            );
            flux_[react->index()][exp].set (GRB_DoubleAttr_Obj, cost);
        }
        void set_biomass_mult (size_t exp, double mult) override
        {
            biomass_[exp].set (GRB_DoubleAttr_Obj, mult);
        }
        
        void finalise_formulation() override
        {
            model_.update();
        }
        
        StatusEnum optimize() override;
        
        SolutionPtr  make_sol (size_t which_experiment) override;
        double get_objective () override
        {
            return model_.get(GRB_DoubleAttr_ObjVal);
        }
        double get_mip_gap () override
        {
            return model_.get(GRB_DoubleAttr_MIPGap);
        }

        void write_model (std::string filename) override
        {
            model_.write (filename);
        }
        size_t num_vars() override
        {
            return (size_t)model_.get (GRB_IntAttr_NumVars);
        }
        size_t num_constraints() override
        {
            return (size_t)model_.get (GRB_IntAttr_NumConstrs);
        }

        // Flux in a single-experiment setting
        GRBVar& flux(size_t k) {return flux_[k][0];}
        GRBVar& flux(size_t k, size_t exp) {return flux_[k][exp];}
        
        void local_set_quiet (bool quiet) override
        {
            model_.set(GRB_IntParam_OutputFlag, quiet ? 0 : 1);
        }
        
        int select_good_reactions (int num_to_select);
        DualValsPtr get_dual_vals();
        double reduced_cost (const Reaction* react);
        
    protected:
        GRBModel model_;

        SolutionPtr cts_sol_;

        // Variable for flux through a reaction 
        std::vector<std::vector<GRBVar>> flux_;
        // Variable for use of a reaction (integer formulation)
        std::vector<GRBVar> use_;
        // The gurobi rep for each constraint
        std::vector<GRBConstr> constr_;

        // The biomass flux var for each experiment
        std::vector<GRBVar> biomass_;
    };
    using GrbLPSolverImpPtr = std::shared_ptr<GrbLPSolverImp>;
}
