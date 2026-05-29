#pragma once

/* An abstraction of the implementation of an LP solver
 */

#include "mosh/gapsolver.h"
#include "mosh/dualvals.h"

namespace mosh
{
    class LPSolverImp 
    {
    public:
        LPSolverImp (Scenario* scenario, Params* params) :
            scenario_(scenario),
            params_(params)
        {
        }

        virtual void init_cts() = 0;
        virtual void init_int(SolutionPtr cts_sol) = 0;

        virtual void make_flux_var (
            const Reaction* react, size_t exp, std::string name,
            double obj_coeff, double lb, double ub
        ) = 0;
        void make_flux_vars (
            const Reaction* react, std::string name,
            double obj_coeff, double lb, double ub
        ) {
            for (size_t exp = 0; exp < num_experiments(); exp++)
                make_flux_var (react, exp, name, obj_coeff, lb, ub);
        }
        virtual void make_use_var (
            const Reaction* react, std::string name, double obj_coeff,
            double lb, double ub
        ) = 0;
        virtual void add_met_constraint (
            const Metabolite* met, size_t experiment, double lb, double ub,
            std::vector<Reaction*>& reacts, std::vector<double>& coeffs
        ) = 0;
        virtual void add_react_link_constraint (
            const Reaction* react, size_t exp
        ) = 0;
        void add_react_link_constraints (const Reaction* react)
        {
            for (size_t exp = 0; exp < num_experiments(); exp++) {
                add_react_link_constraint (react, exp);
            }
        }
        virtual void set_react_bounds (
            const Reaction* react, size_t exp, double lb, double ub
        ) = 0;
        virtual void set_react_cost (
            size_t exp, const Reaction* react, double cost
        ) = 0;
        virtual void set_biomass_mult (size_t exp, double mult) = 0;
        virtual void finalise_formulation() {};
        
        virtual StatusEnum optimize() = 0;

        virtual SolutionPtr make_sol()
        {
            if (num_experiments() > 1) {
                MultiSolPtr sol =
                    std::make_shared<MultiSol> (scenario_, params_);
                for (size_t k = 0; k < num_experiments(); k++)
                    sol->set_sol (k, make_sol(k));
                return sol;
            }
            return make_sol (0);
        }
        
        virtual SolutionPtr make_sol(size_t which_experiment) = 0;
        virtual double get_objective () = 0;
        virtual double get_mip_gap () = 0;
        
        virtual void write_model (std::string filename) = 0;

        virtual size_t num_vars() = 0;
        virtual size_t num_constraints() = 0;

        size_t num_metabolites() const {return scenario_->num_metabolites();}
        const Metabolite* metabolite(size_t k)
        {
            return scenario_->metabolite(k);
        }
        
        size_t num_reactions() const {return scenario_->num_reactions();}
        Reaction* reaction(size_t k) const {return scenario_->reaction(k);}
        
        size_t num_experiments() const {return scenario_->num_experiments();}
        const Experiment* experiment(size_t k) const
        {
            return scenario_->experiment(k);
        }
        
        virtual void local_set_quiet (bool quiet)
        {
        }
        
    protected:
        Scenario* scenario_;
        Params* params_;        
    };
    using LPSolverImpPtr = std::shared_ptr<LPSolverImp>;
}
