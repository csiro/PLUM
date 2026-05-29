#pragma once


#include "lp_lib.h"

#include "mosh/lpsolverimp.h"

#include "lime/numutil.h"

namespace mosh
{
    class LpsLPSolverImp : public LPSolverImp
    {
    public:
        LpsLPSolverImp (Scenario* scenario, Params* params) :
            LPSolverImp (scenario, params),
            model_(NULL),
            add_row_mode_set_(false),
            expect_num_rows_(0),
            lps_coeffs_(new REAL[scenario_->num_reactions()]),
            col_id_(new int[scenario_->num_reactions()]),
            biomass_id_(-1)
        {
        }
        
        virtual ~LpsLPSolverImp()
        {
            DEBUG ('A', "Cleanup");
            if (model_ != NULL)
                delete_lp (model_);
            delete [] lps_coeffs_;
            delete [] col_id_;
        }

        void init_cts() override
        {
            // Set up for continuous prob
            model_ = make_lp(0, scenario_->num_reactions());
            expect_num_rows_ = 2 * scenario_->num_metabolites();
            const char* name = "plum";
            set_lp_name (model_, (char*)name);
            set_verbose (model_, IMPORTANT);
            if (params_->time_limit > 0) {
                long long_limit = (long) params_->time_limit;
                set_timeout (model_, long_limit);
            }
        }

        void init_int(SolutionPtr cts_sol) override
        {
            if (model_ != NULL)
                delete_lp (model_);
            model_ = make_lp(0, 2 * scenario_->num_reactions());
            const char* name = "plum-int";
            set_lp_name (model_, (char*)name);
            set_verbose (model_, IMPORTANT);
            if (params_->time_limit > 0) {
                long long_limit = (long) params_->time_limit;
                set_timeout (model_, long_limit);
            }
            // The number of constraints is
            // - 1 per reaction (linking)
            // - 2 per metabolite per experiment (balance: LB + UB)
            expect_num_rows_ =
                scenario_->num_reactions() + 
                2 * (scenario_->num_metabolites() + scenario_->num_experiments());
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
            return ::get_objective(model_);
        }
        double get_mip_gap () override
        {
            return -1.0f;
        }
        
        void write_model (std::string filename) override
        {
            write_lp (model_, (char*)filename.c_str());
        }

        size_t num_vars() override
        {
            return get_Ncolumns(model_);
        }
        size_t num_constraints() override
        {
            return get_Nrows(model_);
        }

    protected:
        lprec* model_;

        bool add_row_mode_set_;
        int expect_num_rows_;
        REAL *lps_coeffs_;
        int *col_id_;
        
        // Index/ID of biomass flux var for each 
        int biomass_id_;
    };
    using LpsLPSolverImpPtr = std::shared_ptr<LpsLPSolverImp>;
}
