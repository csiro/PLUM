#pragma once

#include "lime/numutil.h"

#include "mosh/lpsolverimp.h"
#include "mosh/highsutil.h"

#include "highs/Highs.h"

namespace mosh
{
    class HighsLPSolverImp : public LPSolverImp
    {
    public:
        HighsLPSolverImp(
            Scenario* scenario, Params* params, int seed
        ) :
            LPSolverImp(scenario, params),
            model_(),
            flux_var_id_(
                scenario_->num_reactions(),
                std::vector<HighsInt> (scenario_->num_experiments())
            ),
            use_var_id_(scenario_->num_reactions()),
            highs_coeffs_(new double[scenario_->num_reactions()]),
            col_id_(new HighsInt[scenario_->num_reactions()]),
            cts_sol_(nullptr),
            biomass_id_(scenario_->num_experiments(), -1)
        {
            if (params->time_limit > 0) {
                check_highs_status (
                    model_.setOptionValue("time_limit", params->time_limit),
                    "setOptionValue-time_limit"
                );
            }
            check_highs_status (
                model_.setOptionValue("threads", params->num_threads),
                "setOptionValue-threads"
            );
            check_highs_status (
                model_.setOptionValue("mip_rel_gap", params->max_mip_gap),
                "setOptionValue-mip_rel_gap"
            );
            check_highs_status (
                model_.setOptionValue(
                    "primal_feasibility_tolerance", params->feasibility_tol
                ),
                "setOptionValue-primal_feasibility_tolerance"
            );
            check_highs_status (
                model_.setOptionValue(
                    "mip_feasibility_tolerance", params->int_feas_tol
                ),
                "setOptionValue-mip_feasibility_tolerance"
            );
            check_highs_status (
                model_.setOptionValue(
                    "random_seed", seed
                ),
                "setOptionValue-random_seed"
            );
        }
        virtual ~HighsLPSolverImp()
        {
            delete [] highs_coeffs_;
            delete [] col_id_;
        }

        void init_cts() override
        {
        }
        
        void init_int(SolutionPtr cts_sol) override
        {
            cts_sol_ = cts_sol;
        }

        void make_flux_var (
            const Reaction* react, size_t exp, std::string name,
            double obj_coeff, double lb, double ub
        ) override;
        void make_use_var (
            const Reaction* react, std::string name, double obj_coeff,
            double lb, double ub
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
        ) override;
        void set_react_cost (
            size_t exp, const Reaction* react, double cost
        ) override;
        void set_biomass_mult (size_t exp, double mult) override;
        
        StatusEnum optimize() override;

        SolutionPtr make_sol(size_t which_experiment) override;
        double get_objective () override
        {
            return model_.getObjectiveValue();
        }
        double get_mip_gap () override
        {
            auto info = model_.getInfo();
            return info.mip_gap;
        }
        
        void write_model (std::string filename) override
        {
            model_.writeModel (filename);
        }

        size_t num_vars() override
        {
            return model_.getNumCol();
        }
        size_t num_constraints() override
        {
            return model_.getNumRow();
        }

        void local_set_quiet (bool quiet) override
        {
            check_highs_status (
                model_.setOptionValue("output_flag", false),
                "setOptionValue-output_flag"
            );
            
        }
        
    protected:
        Highs model_;

        // The flux and use far indices for each reaction in each experiment
        std::vector<std::vector<HighsInt>> flux_var_id_;
        std::vector<HighsInt> use_var_id_;

        /* Working area for calls to "addRow" */
        double *highs_coeffs_;
        HighsInt *col_id_;

        SolutionPtr cts_sol_;
        
        // Index/ID of biomass flux var for each experiment
        std::vector<HighsInt> biomass_id_;
    };
    using HighsLPSolverImpPtr = std::shared_ptr<HighsLPSolverImp>;
}
