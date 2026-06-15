/**
 * @file params.h
 * @brief Algorithm parameters container for metabolic gap-filling optimization
 *
 * This file defines the Params class which holds all configurable parameters
 * used in metabolic network gap-filling algorithms including flux balance analysis,
 * mixed-integer programming, simulated annealing, and tabu search heuristics.
 */
#pragma once

/** A recepticle for all algorithm parameters - calculated and passed
 */

#include <iostream>
#include <set>
#include <memory>

#include "lime/displayable.h"
#include "lime/numutil.h"
#include "lime/config.h"

/**
 * @namespace mosh
 * @brief Namespace for metabolic optimization and simulation heuristics
 */
namespace mosh
{
    // Algorithm parameters. All fields are public
    // This is where to assign defaults (which will be standard across modules)
    
    /** @brief Forward declaration of Scenario class */
    class Scenario;
    /** @brief Forward declaration of Metabolite class */
    class Metabolite;
    /** @brief Forward declaration of Reaction class */
    class Reaction;
  
    /**
     * @class Params
     * @brief Container for all algorithm parameters used in metabolic gap-filling optimization
     *
     * This class serves as a centralized repository for all configurable parameters
     * used across the metabolic optimization algorithms. It includes parameters for:
     * - Mixed-integer programming (MIP) solver settings
     * - Biomass production constraints and search strategies
     * - Metaheuristic algorithms (tabu search, simulated annealing)
     * - Flux thresholds and reaction cost settings
     * - Dummy reaction handling and error detection
     *
     * All member variables are public to allow direct access. Default values are
     * assigned in the constructor and can be overridden via configuration files.
     */
    class Params : public lime::Displayable
    {
    public:
        /**
         * @brief Default constructor initializing all parameters to standard default values
         *
         * Sets default values for integer, floating-point, and boolean parameters used
         * throughout the metabolic optimization algorithms. These defaults are designed
         * to work across different metabolic network models.
         */
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

        /**
         * @brief Finalizes parameter settings after configuration
         *
         * Performs any necessary post-configuration validation or adjustment of parameters.
         * Currently a placeholder for future validation logic.
         */
        void finalise ()
        {
        }
            
        /**
         * @brief Initializes configuration object with current parameter values
         *
         * Populates the provided configuration object with all parameter names and their
         * current values. This is used for configuration file generation and display.
         *
         * @param config Configuration object to populate with parameter values
         */
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

        /**
         * @brief Reads parameter values from configuration object
         *
         * Loads all parameter values from the provided configuration object, overriding
         * default values. Performs validation on certain parameters (e.g., num_intervals >= 1).
         *
         * @param config Configuration object containing parameter values to read
         */
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

        /**
         * @brief Determines if a flux value indicates growth
         *
         * Compares the flux value against the growth_flux_threshold to classify
         * whether the flux represents metabolic growth.
         *
         * @param flux The flux value to evaluate
         * @return true if flux >= growth_flux_threshold, false otherwise
         */
        bool is_growth_flux(double flux) const {
            return flux >= growth_flux_threshold;
        }
        /**
         * @brief Determines if a flux value indicates runaway (unbounded) behavior
         *
         * Compares the flux value against the runaway_flux_threshold to detect
         * unrealistic or unbounded flux values in the metabolic network.
         *
         * @param flux The flux value to evaluate
         * @return true if flux >= runaway_flux_threshold, false otherwise
         */
        bool is_runaway_flux(double flux) const {
            return flux >= runaway_flux_threshold;
        }
        /**
         * @brief Determines if a Biolog score indicates growth
         *
         * Compares the Biolog phenotype score against the growth_biolog_threshold
         * to classify whether the phenotype represents growth.
         *
         * @param biolog The Biolog phenotype score to evaluate
         * @return true if biolog >= growth_biolog_threshold, false otherwise
         */
        bool is_growth_biolog(double biolog) const {
            return biolog >= growth_biolog_threshold;
        }
        /**
         * @brief Checks if a reaction is indicated by genomic evidence
         *
         * Determines whether a reaction has supporting gene annotation based on
         * its cost relative to the gene_ind_cost threshold.
         *
         * @param react Pointer to the reaction to evaluate
         * @return true if the reaction is gene-indicated, false otherwise
         */
        bool is_gene_indicated (const Reaction* react) const;
        /**
         * @brief Checks if a reaction's cost exceeds the maximum allowed cost
         *
         * Determines whether a reaction should be excluded from the initial selection
         * based on its cost exceeding max_react_cost.
         *
         * @param react Pointer to the reaction to evaluate
         * @return true if the reaction cost exceeds max_react_cost, false otherwise
         */
        bool exceeds_max_cost (const Reaction* react) const;
        
        /**
         * @brief Displays all parameter values to an output stream
         *
         * Overrides the Displayable::display method to output all parameter values
         * in a human-readable format via the configuration object.
         *
         * @param os Output stream to write parameter values (default: std::cout)
         */
        virtual void display(std::ostream& os = std::cout) const
        {
            lime::Config config;
            init_config (config);
            config.display (os);
        }


        // Int params
        /** @brief Random number generator seed for reproducibility */
        int seed;
        /** @brief Number of parallel threads for optimization solver */
        int num_threads;
        /** @brief Maximum time limit in seconds for optimization (0 = no limit) */
        int time_limit;
        /** @brief Maximum iterations for biomass objective search algorithm */
        int max_biomass_search_iters;
        /** @brief Maximum iterations for metaheuristic optimization */
        int max_mh_iters;
        /** @brief Extra combinatorial parameter for search space exploration */
        int comb_extra;
        /** @brief Rank constraint type: 0 = none, 1 = relative, 2 = absolute */
        int rank_cons_type;
        // plummx params
        /** @brief Tabu tenure after a reaction becomes incumbent (added or deleted) */
        int incumb_tabu_tenure;
        /** @brief Tabu tenure after a move fails to produce biomass */
        int fail_tabu_tenure;
        /** @brief Tabu tenure after a move does not improve objective */
        int ignore_tabu_tenure;
        /** @brief Number of intervals for target biomass search partitioning */
        int num_intervals;
        /** @brief Number of simulated annealing restarts */
        int sa_restarts;
        /** @brief Maximum iterations without improvement before resetting incumbent to best solution */
        int max_no_improve;

        // Double params
        /** @brief Percentage wiggle room allowed for biomass constraints */
        double biomass_wiggle_pc;
        /** @brief Upper bound on biomass production flux */
        double biomass_ub;
        /** @brief Upper bound on absolute value of objective function */
        double abs_obj_ub;
        /** @brief Initial multiplier for biomass term in objective function */
        double init_biomass_obj_mult;
        /** @brief Interval multiplier for biomass scaling during search */
        double biomass_mult_int;
        /** @brief Optimality multiplier for biomass constraint (biomass >= biomass_opt_mult * max_biomass) */
        double biomass_opt_mult;
        /** @brief Multiplier for increasing biomass objective during initial search phase */
        double biomass_search_mult1;
        /** @brief Multiplier for increasing biomass objective after biomass production is achieved */
        double biomass_search_mult2;
        /** @brief Cost threshold below which a reaction is considered gene-indicated */
        double gene_ind_cost;
        /** @brief Penalty cost assigned to dummy reactions in the gap-filling process */
        double dummy_cost;
        /** @brief Maximum reaction cost threshold; reactions exceeding this are un-selected at start (-1 = ignore) */
        double max_react_cost;
        /** @brief Minimum flux value classified as growth (below this = no growth) */
        double growth_flux_threshold;
        /** @brief Flux value threshold indicating runaway (unbounded) flux behavior */
        double runaway_flux_threshold;
        /** @brief Minimum Biolog phenotype score classified as growth (below this = no growth) */
        double growth_biolog_threshold;
        // Gurobi params
        /** @brief Maximum MIP optimality gap tolerance for Gurobi solver */
        double max_mip_gap;
        /** @brief Feasibility tolerance for constraint satisfaction in Gurobi */
        double feasibility_tol;
        /** @brief Integer feasibility tolerance for integer variables in Gurobi */
        double int_feas_tol;
        /** @brief Small epsilon value for numerical comparisons and tolerances */
        double epsilon;
        // plummx
        /** @brief Multiplier to determine when experimental flux is "over target" */
        double rank_tol_mult;
        /** @brief Multiplier to determine when experimental flux is "excessively" over target */
        double excess_growth_mult;
        //   Simulated Annealing parameters
        /** @brief Target acceptance rate for simulated annealing temperature adjustment */
        double sa_targ_accept;
        /** @brief Target acceptance probability for simulated annealing */
        double sa_targ_prob;
        /** @brief Fraction of iterations with no change triggering simulated annealing termination */
        double sa_frac_no_chg;
        //   Intervals control
        /** @brief Multiplier for calculating lower bound of target interval */
        double lower_target_mult;
        /** @brief Target flux value for interval-based search strategies */
        double target_flux;
        /** @brief Multiplier for tau (relaxation parameter) in PLUMMX objective function */
        double tau_mult;
        /** @brief Multiplier for error term in PLUMMX objective function */
        double error_mult;
        /** @brief Multiplier for reaction cost term in PLUMMX objective function */
        double cost_mult;
        /** @brief Multiplier for runaway flux penalty in PLUMMX objective function */
        double runaways_mult;
        /** @brief Proportion of extra reactions allowed relative to biomass multiplier */
        double mx_biomass_inflation;
        // Adaptation params - signma values
        /** @brief Adaptation sigma parameter for new best solution reward */
        double adapt_sigma1; // New best reward
        /** @brief Adaptation sigma parameter for new incumbent solution */
        double adapt_sigma2; // New incumbant
        /** @brief Adaptation sigma parameter for new Pareto-optimal solution */
        double adapt_sigma3; // New pareto
        /** @brief Learning rate for adaptive parameter adjustment */
        double adapt_learn_rate; // New best reward
        
        // Bool params
        /** @brief Enable absolute value objective function formulation */
        bool use_abs_obj;
        /** @brief Enable error term in objective function */
        bool use_error;
        /** @brief Enable dummy reactions for gap-filling */
        bool use_dummy;
        /** @brief Enable dummy biomass reaction for demand metabolites */
        bool use_dummy_biomass_react_dm;
        /** @brief Enable dummy biomass reaction for exchange reactions */
        bool use_dummy_biomass_react_ex;
        /** @brief Enable dummy biomass production for demand metabolites */
        bool use_dummy_biomass_prod_dm;
        /** @brief Preserve dummy reactions in final solution (don't remove) */
        bool preserve_dummies;
        /** @brief Enable detection and penalization of runaway fluxes */
        bool detect_runaway;
        /** @brief Use unit cost (1.0) for all reactions instead of database costs */
        bool unit_cost;
    };
    
    /** @brief Shared pointer type alias for Params objects */
    using ParamsPtr = std::shared_ptr<Params>;
}
