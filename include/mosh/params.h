#pragma once

/** A recepticle for all algorithm parameters - calculated and passed
 */

#include <iostream>
#include <set>
#include <memory>

#include "lime/displayable.h"
#include "lime/numutil.h"
#include "lime/config.h"

namespace mosh
{
    // Algorithm parameters. All fields are public
    // This is where to assign defaults (which will be standard across modules)
    
    class Scenario;
    class Metabolite;
    class Reaction;
  
    class Params : public lime::Displayable
    {
    public:
        Params () :
            seed(0),
            num_threads(8),
            time_limit(0),
            max_biomass_search_iters(20),
            max_mh_iters(100),
            comb_extra(0),
            rank_cons_type(1),
            incumb_tabu_tenure(500),
            fail_tabu_tenure(200),
            ignore_tabu_tenure(25),
            num_intervals(1),
            sa_restarts(0),
            max_no_improve(1000),
            
            biomass_wiggle_pc(0.0f),
            biomass_ub(0.0f),
            abs_obj_ub(0.0f),
            init_biomass_obj_mult(100.0f),
            biomass_mult_int(10.0f),
            biomass_opt_mult(0.9999f),
            biomass_search_mult1(2.0f),
            biomass_search_mult2(1.25f),
            gene_ind_cost(1.0f),
            dummy_cost(1.0e4f),
            max_react_cost(-1.0f),
            growth_flux_threshold (0.01f),
            runaway_flux_threshold(10.0f),
            growth_biolog_threshold(100.0f),
            max_mip_gap(0.01f),
            feasibility_tol(1e-8f),
            int_feas_tol(1e-8f),
            epsilon(1e-3f),
            rank_tol_mult(1.01f),
            excess_growth_mult(1.5f),
            sa_targ_accept(1.1f),
            sa_targ_prob(0.25f),
            sa_frac_no_chg(0.05),
            lower_target_mult(0.5f),
            target_flux(0.0f),
            tau_mult(10.0f),
            error_mult(1.0f), 
            cost_mult(10.0f), 
            runaways_mult(10000.0f), 
            mx_biomass_inflation(1.25f), 
            adapt_sigma1(10.0f), 
            adapt_sigma2(3.0f), 
            adapt_sigma3(1.0f), 
            adapt_learn_rate(0.3f), 
            
            use_abs_obj(false),
            use_error(false),
            use_dummy(false),
            use_dummy_biomass_react_dm(false),
            use_dummy_biomass_react_ex(false),
            use_dummy_biomass_prod_dm(false),
            preserve_dummies(false),
            detect_runaway(true),
            unit_cost(false)
        {
        }

        void finalise ()
        {
        }
            
        void init_config (lime::Config& config) const
        {
            // Int config
            config.addItem ("seed", seed);
            config.addItem ("threads", num_threads);
            config.addItem ("timelimit", time_limit);
            config.addItem ("max_biomass_iters", max_biomass_search_iters);
            config.addItem ("max_mh_iters", max_mh_iters);
            config.addItem ("comb_extra", comb_extra);
            config.addItem ("rank_cons_type", rank_cons_type);
            config.addItem ("incumb_tabu_tenure", incumb_tabu_tenure);
            config.addItem ("fail_tabu_tenure", fail_tabu_tenure);
            config.addItem ("ignore_tabu_tenure", ignore_tabu_tenure);
            config.addItem ("num_intervals", num_intervals);
            config.addItem ("sa_restarts", sa_restarts);
            config.addItem ("max_no_improve", max_no_improve);

            // Double config
            config.addItem ("biomass_wiggle_pc", biomass_wiggle_pc);
            config.addItem ("biomass_ub", biomass_ub);
            config.addItem ("abs_obj_ub", abs_obj_ub);
            config.addItem ("init_biomass_mult", init_biomass_obj_mult);
            config.addItem ("biomass_mult_int", biomass_mult_int);
            config.addItem ("biomass_opt_mult", biomass_opt_mult);
            config.addItem ("biomass_search_mult1", biomass_search_mult1);
            config.addItem ("biomass_search_mult2", biomass_search_mult2);
            config.addItem ("gene_ind_cost", gene_ind_cost);
            config.addItem ("dummy_cost", dummy_cost);
            config.addItem ("max_react_cost", max_react_cost);
            config.addItem ("growth_biolog_threshold", growth_biolog_threshold);
            config.addItem ("max_mip_gap", max_mip_gap);
            config.addItem ("epsilon", epsilon);
            config.addItem ("rank_tol_mult", rank_tol_mult);
            config.addItem ("excess_growth_mult", excess_growth_mult);
            config.addItem ("sa_targ_accept", sa_targ_accept);
            config.addItem ("sa_targ_prob", sa_targ_prob);
            config.addItem ("sa_frac_no_chg", sa_frac_no_chg);
            config.addItem ("lower_target_mult", lower_target_mult);
            config.addItem ("target_flux", target_flux);
            config.addItem ("tau_mult", tau_mult);
            config.addItem ("error_mult", error_mult);
            config.addItem ("cost_mult", cost_mult);
            config.addItem ("runaways_mult", runaways_mult);
            config.addItem ("mx_biomass_inflation", mx_biomass_inflation);
            config.addItem ("adapt_sigma1", adapt_sigma1);
            config.addItem ("adapt_sigma2", adapt_sigma2);
            config.addItem ("adapt_sigma3", adapt_sigma3);
            config.addItem ("adapt_learn_rate", adapt_learn_rate);
            
            // Bool config
            config.addItem ("use_abs", use_abs_obj);
            config.addItem ("use_error", use_error);
            config.addItem ("use_dummy", use_dummy);
            config.addItem (
                "use_dummy_biomass_react_dm", use_dummy_biomass_react_dm
            );
            config.addItem (
                "use_dummy_biomass_react_ex", use_dummy_biomass_react_ex
            );
            config.addItem (
                "use_dummy_biomass_prod_dm", use_dummy_biomass_prod_dm
            );
            config.addItem ("preserve_dummies", preserve_dummies);
            config.addItem ("detect_runaway", detect_runaway);
            config.addItem ("unit_cost", unit_cost);
        }

        void read_config (lime::Config& config) 
        {
            // Int config
            seed = config.getInt ("seed", seed);
            num_threads = config.getInt ("threads", num_threads);
            time_limit = config.getInt ("timelimit", time_limit);
            max_biomass_search_iters =
                config.getInt ("max_biomass_iters", max_biomass_search_iters);
            max_mh_iters = config.getInt ("max_mh_iters", max_mh_iters);
            comb_extra = config.getInt ("comb_extra", comb_extra);
            rank_cons_type = config.getInt ("rank_cons_type", rank_cons_type);
            incumb_tabu_tenure =
                config.getInt ("incumb_tabu_tenure", incumb_tabu_tenure);
            fail_tabu_tenure =
                config.getInt ("fail_tabu_tenure", fail_tabu_tenure);
            ignore_tabu_tenure =
                config.getInt ("ignore_tabu_tenure", ignore_tabu_tenure);
            num_intervals =
                config.getInt ("num_intervals", num_intervals);
            if (num_intervals < 1)
                num_intervals = 1;
            sa_restarts = config.getInt ("sa_restarts", sa_restarts);
            max_no_improve = config.getInt ("max_no_improve", max_no_improve);

            // Double config
            biomass_wiggle_pc =
                config.getDouble ("biomass_wiggle_pc", biomass_wiggle_pc);
            biomass_ub =
                config.getDouble ("biomass_ub", biomass_ub);
            abs_obj_ub =
                config.getDouble ("abs_obj_ub", abs_obj_ub);
            init_biomass_obj_mult =
                config.getDouble ("init_biomass_mult", init_biomass_obj_mult);
            biomass_mult_int =
                config.getDouble ("biomass_mult_int", biomass_mult_int);
            biomass_opt_mult =
                config.getDouble ("biomass_opt_mult", biomass_opt_mult);
            biomass_search_mult1 =
                config.getDouble ("biomass_search_mult1", biomass_search_mult1);
            biomass_search_mult2 =
                config.getDouble ("biomass_search_mult2", biomass_search_mult2);
            gene_ind_cost =
                config.getDouble ("gene_ind_cost", gene_ind_cost);
            dummy_cost =
                config.getDouble ("dummy_cost", dummy_cost);
            max_react_cost =
                config.getDouble ("max_react_cost", max_react_cost);
            growth_biolog_threshold =
                config.getDouble (
                    "growth_biolog_threshold", growth_biolog_threshold
                );
            max_mip_gap = config.getDouble ("max_mip_gap", max_mip_gap);
            epsilon = config.getDouble ("epsilon", epsilon);
            rank_tol_mult = config.getDouble ("rank_tol_mult", rank_tol_mult);
            excess_growth_mult =
                config.getDouble ("excess_growth_mult", excess_growth_mult);
            sa_targ_accept =
                config.getDouble ("sa_targ_accept", sa_targ_accept);
            sa_targ_prob = config.getDouble ("sa_targ_prob", sa_targ_prob);
            sa_frac_no_chg =
                config.getDouble ("sa_frac_no_chg", sa_frac_no_chg);
            lower_target_mult =
                config.getDouble ("lower_target_mult", lower_target_mult);
            target_flux =
                config.getDouble ("target_flux", target_flux);
            tau_mult =
                config.getDouble ("tau_mult", tau_mult);
            error_mult =
                config.getDouble ("error_mult", error_mult);
            cost_mult =
                config.getDouble ("cost_mult", cost_mult);
            runaways_mult =
                config.getDouble ("runaways_mult", runaways_mult);
            mx_biomass_inflation =
                config.getDouble ("mx_biomass_inflation", mx_biomass_inflation);
            adapt_sigma1 =
                config.getDouble ("adapt_sigma1", adapt_sigma1);
            adapt_sigma2 =
                config.getDouble ("adapt_sigma2", adapt_sigma2);
            adapt_sigma3 =
                config.getDouble ("adapt_sigma3", adapt_sigma3);
            adapt_learn_rate =
                config.getDouble ("adapt_learn_rate", adapt_learn_rate);
            
            // Bool config
            use_abs_obj = config.getBool ("use_abs", use_abs_obj);
            use_error = config.getBool ("use_error", use_error);
            use_dummy = config.getBool ("use_dummy", use_dummy);
            use_dummy_biomass_react_dm =
                config.getBool (
                    "use_dummy_biomass_react_dm", use_dummy_biomass_react_dm
                );
            use_dummy_biomass_react_ex =
                config.getBool (
                    "use_dummy_biomass_react_ex", use_dummy_biomass_react_ex
                );
            use_dummy_biomass_prod_dm =
                config.getBool (
                    "use_dummy_biomass_prod_dm", use_dummy_biomass_prod_dm
                );
            preserve_dummies =
                config.getBool ("preserve_dummies", preserve_dummies);
            detect_runaway = config.getBool ("detect_runaway", detect_runaway);
            unit_cost = config.getBool ("unit_cost", unit_cost);
        }

        bool is_growth_flux(double flux) const {
            return flux >= growth_flux_threshold;
        }
        bool is_runaway_flux(double flux) const {
            return flux >= runaway_flux_threshold;
        }
        bool is_growth_biolog(double biolog) const {
            return biolog >= growth_biolog_threshold;
        }
        bool is_gene_indicated (const Reaction* react) const;
        bool exceeds_max_cost (const Reaction* react) const;
        
        virtual void display(std::ostream& os = std::cout) const
        {
            lime::Config config;
            init_config (config);
            config.display (os);
        }


        // Int params
        int seed;
        int num_threads;
        int time_limit;
        int max_biomass_search_iters;
        int max_mh_iters;
        int comb_extra;
        // Rank cons type 0 = none, 1 = relative, 2 = absolute
        int rank_cons_type;
        // plummx params
        //   Tabu tenure after becoming incumbant (add or delete)
        int incumb_tabu_tenure;
        //   Tabu tenure after failing (no biomass)
        int fail_tabu_tenure;
        //   Tabu tenure after not improving
        int ignore_tabu_tenure;
        //   Number of intervals looking for a target biomass
        int num_intervals;
        int sa_restarts;
        // Numer of iters in plummx with no improve (then reset incumb to best)
        int max_no_improve;

        // Double params
        double biomass_wiggle_pc;
        double biomass_ub;
        double abs_obj_ub;
        double init_biomass_obj_mult;
        double biomass_mult_int;
        // Optimility of biomass for formulation 3
        // (Biomass is constrained to be >= biomass_opt_mult * max biomass)
        double biomass_opt_mult;
        // Amount (multiple) biomass mult is increased in biomass search
        double biomass_search_mult1;
        // Amount (multiple) biomass mult is increased in biomass search after
        // biomass has been produced
        double biomass_search_mult2;
        // Cost that implies the reaction is "gene-inidicated"
        double gene_ind_cost;
        // Cost used for dummy reactions
        double dummy_cost;
        // Max react cost - higher is un-selected at start (-1 = ignore)
        double max_react_cost;
        // What is the max flux still classified as "no growth"
        double growth_flux_threshold;
        // What is the flux classified as "runaway"
        double runaway_flux_threshold;
        // What is the max biolog score still classified as "no growth"
        double growth_biolog_threshold;
        // Gurobi params
        double max_mip_gap;
        double feasibility_tol;
        double int_feas_tol;
        double epsilon;
        // plummx 
        //   Multiplier to indicate when an experiment is "over target"
        double rank_tol_mult;
        //   Multiplier to indicate when an experiment is "excessively" over
        //   target
        double excess_growth_mult;
        //   Simulated Annealing parameters
        double sa_targ_accept;
        double sa_targ_prob;
        double sa_frac_no_chg;
        //   Intervals control
        double lower_target_mult;
        double target_flux;
        // Multiplier for tau in objectives for plummx
        double tau_mult;
        // Multiplier for error in objectives for plummx
        double error_mult;
        // Multiplier for cost in objectives for plummx
        double cost_mult;
        // Multiplier for runaways in objectives for plummx
        double runaways_mult;
        // Multiplier for biomass mult - proportion of extra reactions allowed
        double mx_biomass_inflation;
        // Adaptation params - signma values
        double adapt_sigma1; // New best reward
        double adapt_sigma2; // New incumbant
        double adapt_sigma3; // New pareto
        double adapt_learn_rate; // New best reward
        
        // Bool params
        bool use_abs_obj;
        bool use_error;
        bool use_dummy;
        bool use_dummy_biomass_react_dm;
        bool use_dummy_biomass_react_ex;
        bool use_dummy_biomass_prod_dm;
        bool preserve_dummies;
        bool detect_runaway;
        bool unit_cost;
    };
    
    using ParamsPtr = std::shared_ptr<Params>;
}
