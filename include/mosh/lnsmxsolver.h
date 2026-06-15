/**
 * @file lnsmxsolver.h
 * @brief Large Neighborhood Search Mixed solver for metabolic gap-filling
 *
 * This file defines the LnsMxSolver class, which implements a Large Neighborhood Search (LNS)
 * metaheuristic with mixed objectives for solving metabolic gap-filling problems. The solver
 * uses simulated annealing and Pareto optimization to find optimal sets of reactions that
 * minimize cost, error, and Kendall's tau while satisfying flux balance constraints.
 */
#pragma once

#include <vector>
#include <memory>
#include <sstream>

#include "lime/rand.h"
#include "lime/debug.h"
#include "lime/biaschoicet.h"
#include "lime/simanneal.h"
#include "lime/adaptchoice.h"
#include "lime/pareto.h"

#include "mosh/gapsolver.h"
#include "mosh/lpsolver.h"
#include "mosh/lpsolverimp.h"
#ifdef PLUM_GUROBI
#include "mosh/grblpsolverimp.h"
#endif
#include "mosh/multisol.h"

/**
 * @brief Namespace for metabolic optimization and solver hierarchy
 */
namespace mosh
{
    /**
     * @class LnsMxSolver
     * @brief Large Neighborhood Search solver with mixed multi-objective optimization
     *
     * Implements a metaheuristic approach to metabolic gap-filling that uses Large Neighborhood Search
     * combined with simulated annealing and Pareto optimization. The solver iteratively adds/removes
     * reactions and adjusts costs to find solutions that balance multiple objectives: minimizing reaction
     * cost, minimizing flux error against target data, and maximizing Kendall's tau correlation.
     *
     * The solver supports various objective function combinations and uses adaptive selection of
     * neighborhood search methods during optimization.
     */
    class LnsMxSolver : public GapSolver
    {
    public:

        /** @brief Type alias for list of reaction pointers */
        using ReactList = std::list<Reaction*>;
        /** @brief Type alias for vector of boolean flags for reactions */
        using ReactBools = std::vector<bool>;
        
        /**
         * @enum ObjElt
         * @brief Elements that can compose a multi-objective optimization
         */
        enum ObjElt {COST, TAU, ERROR};
        /** @brief Names of objective elements for display */
        std::vector<const char*> objelt_name_ = {
            "COST", "TAU", "ERROR"
        };
        
        /**
         * @enum WhichObj
         * @brief Enumeration of available objective function combinations
         *
         * Defines different combinations of cost, tau, and error objectives.
         * The naming convention uses '__' to separate primary and secondary objectives.
         */
        enum WhichObj {
            ERROR__COST = 0,
            TAU_ERROR__COST = 1,
            TAU_COST__ERROR = 2,
            COST_TAU__ERROR = 3,
            ERROR_COST__TAU = 4,
            ALL_OBJ = 5,
            TAU__COST = 6,
            COST__TAU = 7,
            NUM_OBJ
        };
        /** @brief Names of objective combinations for display */
        std::vector<const char*> obj_name_ = {
            "ERROR/COST", "TAU+ERROR/COST", "TAU+COST/ERROR", "COST+TAU/ERROR",
            "ERROR+COST/TAU", "TAU+ERROR+COST/COST", "TAU/COST", "COST/TAU"
        };
        
        /**
         * @class SolAttributes
         * @brief Container for computed attributes and metrics of a solution
         *
         * Stores and computes various quality metrics for a multi-solution including
         * objective values, error metrics, cost information, growth predictions, and
         * reaction usage statistics. Used for solution comparison and Pareto optimization.
         */
        class SolAttributes
        {
        public:
            /**
             * @brief Construct solution attributes from a multi-solution
             * @param lnsmx Pointer to the parent LnsMxSolver
             * @param sol Multi-solution to compute attributes for
             * @param non_unit_react Boolean vector indicating which reactions have non-unit cost
             */
            SolAttributes (
                LnsMxSolver* lnsmx, MultiSolPtr sol,
                ReactBools& non_unit_react
            ) :
                sol_(sol),
                obj_(0.0f),
                obj2_(0.0f),
                multi_obj_({0.0f,0.0f,0.0f}),
                growth_mismatch_(0),
                pc_true_pos_(0),
                pc_true_neg_(0),
                error_(lnsmx->mean_squared_error (sol)),
                tau_(lnsmx->calc_kendall_tau (sol)),
                cost_(sol->abs_obj_value()),
                cost_ratio_(cost_ / lnsmx->base_cost_),
                max_cost_(sol->max_cost()),
                num_runaways_(sol->num_runaways()),
                rank0_target_(
                    lnsmx->target_flux(lnsmx->scenario_->biolog_rank0())
                ),
                unselected_(),
                non_unit_react_(non_unit_react),
                num_non_unit_(0),
                non_unit_used_(0)
            {
                calc_growth_mismatch (lnsmx);
                calc_primary_obj (lnsmx);
                calc_secondary_obj (lnsmx);
                calc_multi_obj (lnsmx);
                for (auto& react : lnsmx->scenario_->reactions()) {
                    if (!react->is_selected())
                        unselected_.push_back (react.get());
                }
                for (size_t k = 0; k < non_unit_react_.size(); k++) {
                    if (non_unit_react_[k]) {
                        num_non_unit_++;
                        if (sol->uses_react(k))
                            non_unit_used_++;
                    }
                }
            }

            /**
             * @brief Copy constructor
             * @param other Source SolAttributes object to copy from
             */
            SolAttributes (SolAttributes& other) :
                sol_(other.sol_),
                obj_(other.obj_),
                obj2_(other.obj2_),
                multi_obj_(other.multi_obj_),
                growth_mismatch_(other.growth_mismatch_),
                pc_true_pos_(other.pc_true_pos_),
                pc_true_neg_(other.pc_true_neg_),
                error_(other.error_),
                tau_(other.tau_),
                cost_(other.cost_),
                cost_ratio_(other.cost_ratio_),
                max_cost_(other.max_cost_),
                num_runaways_(other.num_runaways_),
                rank0_target_(other.rank0_target_),
                unselected_(other.unselected_.begin(), other.unselected_.end()),
                non_unit_react_(other.non_unit_react_),
                num_non_unit_(other.num_non_unit_),
                non_unit_used_(other.non_unit_used_)
            {
            }

            /**
             * @brief Get the multi-solution pointer
             * @return Reference to the MultiSolPtr
             */
            MultiSolPtr& sol() {return sol_;}
            /**
             * @brief Get individual solution by index
             * @param k Index of the solution
             * @return Pointer to the k-th solution
             */
            SolutionPtr sol(size_t k) {return sol_->sol(k);}
            /**
             * @brief Get primary objective value
             * @return Primary objective value
             */
            double obj() const {return obj_;}
            /**
             * @brief Get secondary objective value
             * @return Secondary objective value
             */
            double obj2() const {return obj_;}
            /**
             * @brief Get multi-objective vector for Pareto optimization
             * @return Reference to vector of objective values
             */
            std::vector<double>& multi_obj() {return multi_obj_;}
            /**
             * @brief Get count of growth prediction mismatches
             * @return Number of experiments where predicted growth disagrees with expected
             */
            int growth_mismatch() const {return growth_mismatch_;}
            /**
             * @brief Get percentage of true positive growth predictions
             * @return Percentage of growth experiments correctly predicted
             */
            int pc_true_positive() const {return pc_true_pos_;}
            /**
             * @brief Get percentage of true negative growth predictions
             * @return Percentage of no-growth experiments correctly predicted
             */
            int pc_true_negative() const {return pc_true_neg_;}
            /**
             * @brief Get mean squared error against target fluxes
             * @return Mean squared error value
             */
            double error() const {return error_;}
            /**
             * @brief Get Kendall's tau correlation coefficient
             * @return Kendall's tau value
             */
            double tau() const {return tau_;}
            /**
             * @brief Get total absolute cost of selected reactions
             * @return Total reaction cost
             */
            double cost() const {return cost_;}
            /**
             * @brief Get cost ratio relative to base cost
             * @return Cost divided by base cost
             */
            double cost_ratio() const {return cost_ratio_;}
            /**
             * @brief Get maximum individual reaction cost
             * @return Maximum cost among selected reactions
             */
            double max_cost() const {return max_cost_;}
            /**
             * @brief Get count of runaway reactions (unbounded fluxes)
             * @return Number of runaway reactions
             */
            int num_runaways() const {return num_runaways_;}
            /**
             * @brief Get target flux for rank 0 experiment
             * @return Target flux value
             */
            double rank0_target() const {return rank0_target_;}
            /**
             * @brief Get list of unselected reactions
             * @return Reference to list of unselected reaction pointers
             */
            ReactList& unselected() {return unselected_;}
            /**
             * @brief Get count of unselected reactions
             * @return Number of unselected reactions
             */
            size_t num_unselected() {return unselected_.size();}
            /**
             * @brief Get boolean vector of non-unit cost reactions
             * @return Reference to boolean vector
             */
            ReactBools& non_unit_react() {return non_unit_react_;}
            /**
             * @brief Get count of reactions with non-unit cost
             * @return Number of non-unit cost reactions
             */
            size_t num_non_unit() const {return num_non_unit_;}
            /**
             * @brief Get count of non-unit cost reactions used in solution
             * @return Number of non-unit cost reactions with non-zero flux
             */
            size_t non_unit_used() const {return non_unit_used_;}

            /**
             * @brief Get biomass flux from solution
             * @return Biomass flux value, or 0 if solution is null
             */
            double biomass_flux () const
            {
                return sol_ == nullptr ? 0.0f : sol_->biomass_flux();
            }

            /**
             * @brief Calculate primary objective value based on solver configuration
             * @param lnsmx Pointer to the parent solver to access objective settings
             */
            void calc_primary_obj (const LnsMxSolver* lnsmx)
            {
                obj_ = 0.0f;
                switch (lnsmx->which_obj_) {
                case ERROR__COST:
                    // Primary is only error
                    obj_ = error_;
                    break;
                case TAU_ERROR__COST:
                    // Primary is Tau then error
                    obj_ =
                        lnsmx->params_->tau_mult * (1.0f - tau_) +
                        error_;
                    break;
                case TAU_COST__ERROR:
                    // Primary is Tau then cost
                    obj_ =
                        lnsmx->params_->tau_mult * (1.0f - tau_) +
                        cost_ratio_;
                    break;
                case COST_TAU__ERROR:
                    // Primary is cost then Tau
                    obj_ =
                        lnsmx->params_->cost_mult * cost_ratio_ +
                        (1.0f - tau_);
                    break;
                case ERROR_COST__TAU:
                    // Primary is cost then Tau
                    obj_ =
                        lnsmx->params_->error_mult * error_ +
                        cost_ratio_;
                    break;
                case ALL_OBJ:
                    // Primary is Tau then error, then cost
                    obj_ =
                        lnsmx->params_->tau_mult * (1.0f - tau_) +
                        lnsmx->params_->error_mult * error_ +
                        lnsmx->params_->cost_mult * cost_ratio_;
                    break;
                case TAU__COST:
                    // Primary is tau
                    obj_ = (1.0f - tau_);
                    break;
                case COST__TAU:
                    // Primary is cost
                    obj_ = cost_ratio_;
                    break;
                }
                // Whichever obj, also add runaways cost (should dominate)
                obj_ += num_runaways_ * lnsmx->params_->runaways_mult;
                DEBUG ('x', "Primary obj " << obj_);
            }

            /**
             * @brief Calculate secondary objective value based on solver configuration
             * @param lnsmx Pointer to the parent solver to access objective settings
             */
            void calc_secondary_obj (const LnsMxSolver* lnsmx)
            {
                obj2_ = 0.0f;
                switch (lnsmx->which_obj_) {
                case ERROR__COST:
                case TAU_ERROR__COST:
                case ALL_OBJ:
                case TAU__COST:
                    // Second is cost
                    obj2_ = cost_;
                    break;
                case TAU_COST__ERROR:
                case COST_TAU__ERROR:
                    // Second is error
                    obj2_ = error_;
                    break;
                case COST__TAU:
                case ERROR_COST__TAU:
                    // Second is tau (maximised)
                    obj2_ = 1.0f - tau_;
                    break;
                }
                DEBUG ('x', "Secondary obj " << obj2_);
            }

            /**
             * @brief Calculate multi-objective vector for Pareto front
             * @param lnsmx Pointer to the parent solver to access objective configuration
             */
            void calc_multi_obj (const LnsMxSolver* lnsmx)
            {
                for (size_t k = 0; k < multi_obj_.size(); k++) {
                    switch (lnsmx->obj_elt(k)) {
                    case COST:   multi_obj_[k] = cost(); break;
                    case TAU:   multi_obj_[k] = 1.0f - tau();   break;
                    case ERROR: multi_obj_[k] = error(); break;
                    }
                }
            }

            /**
             * @brief Calculate growth prediction accuracy metrics
             *
             * Computes the number of mismatches between predicted and expected growth,
             * as well as true positive and true negative percentages.
             *
             * @param lnsmx Pointer to the parent solver to access experiment data
             */
            void calc_growth_mismatch (const LnsMxSolver* lnsmx)
            {
                // Count number of cases where required
                // 'growth' media does not grow.
                // We don't care about no-growth media - it may or may not
                growth_mismatch_ = 0;
                int pos_match = 0;
                int num_growth = 0;
                int neg_match = 0;
                int num_nogrowth = 0;

                for (int k = 0; k < sol_->num_exp(); k++) {
                    DEBUG (
                        'x', "  Sol " << k <<
                        " " << lnsmx->experiment(k)->name() <<
                        " reqd growth " << lnsmx->experiment(k)->is_growth() <<
                        " is growth " <<
                        lnsmx->params_->is_growth_flux(sol_->biomass_flux(k))
                    );
                    bool is_growth_flux =
                        lnsmx->params_->is_growth_flux(sol_->biomass_flux(k));

                    if (lnsmx->experiment(k)->is_growth()) {
                        num_growth++;
                        if (is_growth_flux)
                            pos_match++;
                        else
                            growth_mismatch_++;
                    }
                    else {
                        num_nogrowth++;
                        if (!is_growth_flux)
                            neg_match++;
                    }
                }
                pc_true_pos_ = 100;
                if (num_growth > 0)
                    pc_true_pos_ = (int) round(100.0f * pos_match / num_growth);
                pc_true_neg_ = 100;
                if (num_nogrowth > 0)
                    pc_true_neg_ =
                        (int) round(100.0f * neg_match / num_nogrowth);
            }

            /**
             * @brief Convert solution attributes to string representation
             * @return String containing all solution metrics
             */
            std::string to_string() {
                std::stringstream str;
                str <<
                    "obj " << obj_ <<
                    " obj2 " << obj2_ <<
                    " multi_obj { " << multi_obj_[0] <<
                    " " << multi_obj_[1] <<
                    " " << multi_obj_[2] << " }" <<
                    " error " << error_ <<
                    " tau " << tau_ <<
                    " totcost " << cost_ <<
                    " cost_ratio " << cost_ratio_ <<
                    " max_cost " << max_cost_ <<
                    " runaways " << num_runaways_ <<
                    " growth_mismatch " << growth_mismatch_ <<
                    " truepospc " << pc_true_pos_ <<
                    " truenegpc " << pc_true_neg_ <<
                    " unsel " << unselected_.size() <<
                    " num_nonunit " << num_non_unit_ <<
                    " nonunit_used " << non_unit_used_;
                return str.str();
            }
            
        private:
            MultiSolPtr sol_; /**< Pointer to the multi-solution */
            double obj_; /**< Primary objective value */
            double obj2_; /**< Secondary objective value */
            std::vector<double> multi_obj_; /**< Multi-objective vector for Pareto optimization */
            int growth_mismatch_; /**< Number of growth prediction mismatches */
            int pc_true_pos_; /**< Percentage of true positive growth predictions */
            int pc_true_neg_; /**< Percentage of true negative growth predictions */
            double error_; /**< Mean squared error against target fluxes */
            double tau_; /**< Kendall's tau correlation coefficient */
            double cost_; /**< Total absolute cost of reactions */
            double cost_ratio_; /**< Cost ratio relative to base cost */
            double max_cost_; /**< Maximum individual reaction cost */
            int num_runaways_; /**< Count of runaway reactions */
            double rank0_target_; /**< Target flux for rank 0 experiment */
            ReactList unselected_; /**< List of unselected reactions */
            ReactBools non_unit_react_; /**< Boolean flags for non-unit cost reactions */
            size_t num_non_unit_; /**< Count of non-unit cost reactions */
            size_t non_unit_used_; /**< Count of non-unit cost reactions used */
        };
        /** @brief Shared pointer type for SolAttributes */
        using SolAttributesPtr = std::shared_ptr<SolAttributes>;        
        
        /**
         * @brief Construct a Large Neighborhood Search solver
         * @param scenario Pointer to the metabolic scenario
         * @param params Pointer to solver parameters
         * @param flavour Solution flavor (e.g., minimize/maximize biomass)
         * @param which_obj Which objective function combination to use
         * @param seed Random number generator seed
         * @param max_iters Maximum number of iterations
         * @param make_all_unit Whether to force all reactions to unit cost
         * @param progname Program name for logging
         */
        LnsMxSolver(
            Scenario* scenario, Params* params,
            Flavour flavour, WhichObj which_obj,
            int seed, int max_iters, bool make_all_unit, std::string progname
        );

        /**
         * @brief Main solve method implementing LNS algorithm
         * @return Pointer to best solution found
         */
        virtual SolutionPtr solve() override;
        /**
         * @brief Solve by iteratively adding reactions
         * @return Pointer to solution
         */
        SolutionPtr solve_by_add();

        /**
         * @brief Get which objective function combination is being used
         * @return WhichObj enumeration value
         */
        WhichObj which_obj() const {return which_obj_;}
        /**
         * @brief Get objective element at index k
         * @param k Index into objective elements
         * @return ObjElt enumeration value
         */
        ObjElt obj_elt (size_t k) const {return obj_elt_[k];}
        
        /**
         * @brief Get the best multi-solution found
         * @return Pointer to best multi-solution, or nullptr if none
         */
        MultiSolPtr best_sol() const
        {
            return best_sol_ == nullptr ? nullptr : best_sol_->sol();
        }
        /**
         * @brief Get attributes of the best solution
         * @return Shared pointer to best solution attributes
         */
        SolAttributesPtr best_attr() {return best_sol_;}

        /**
         * @brief Calculate Kendall's tau correlation coefficient for a solution
         * @param sol Multi-solution to evaluate
         * @return Kendall's tau value
         */
        double calc_kendall_tau (MultiSolPtr sol) const;
        /**
         * @brief Calculate mean squared error against target fluxes
         * @param sol Multi-solution to evaluate
         * @return Mean squared error value
         */
        double mean_squared_error (MultiSolPtr sol) const;
        /**
         * @brief Get number of reactions used in best solution
         * @return Count of used reactions
         */
        int num_used () const;
        /**
         * @brief Get number of non-unit cost reactions
         * @return Count of non-unit cost reactions
         */
        size_t num_non_unit() const
        {
            return best_sol_ == nullptr ? 0 : best_sol_->num_non_unit();
        }
        /**
         * @brief Get number of non-unit cost reactions used in best solution
         * @return Count of non-unit cost reactions with non-zero flux
         */
        size_t non_unit_used() const
        {
            return best_sol_ == nullptr ? 0 : best_sol_->non_unit_used();
        }
        /**
         * @brief Get percentage of true positive growth predictions from best solution
         * @return Percentage (0-100) of correctly predicted growth experiments
         */
        int pc_true_positive() const
        {
            return best_sol_ == nullptr ? 0 : best_sol_->pc_true_positive();
        }
        /**
         * @brief Get percentage of true negative growth predictions from best solution
         * @return Percentage (0-100) of correctly predicted no-growth experiments
         */
        int pc_true_negative() const
        {
            return best_sol_ == nullptr ? 0 : best_sol_->pc_true_negative();
        }
        
        /**
         * @brief Write list of unselected reactions to output stream
         * @param out Output stream
         */
        void write_best_unselected (std::ostream& out);
        /**
         * @brief Write best solution to output stream
         * @param out Output stream
         */
        void write_best_sol (std::ostream& out);
        
        /**
         * @brief Get filename for saving best solution
         * @return Filename string
         */
        std::string save_best_fn() const {return save_best_fn_;}
        /**
         * @brief Set filename for saving best solution
         * @param save_best_fn Filename to use
         */
        void set_save_best_fn (std::string save_best_fn) {
            save_best_fn_ = save_best_fn;
        }
        /**
         * @brief Get filename for saving best solution details
         * @return Filename string
         */
        std::string save_best_sol_fn() const {return save_best_sol_fn_;}
        /**
         * @brief Set filename for saving best solution details
         * @param save_best_sol_fn Filename to use
         */
        void set_save_best_sol_fn (std::string save_best_sol_fn) {
            save_best_sol_fn_ = save_best_sol_fn;
        }
        /**
         * @brief Set filename for progress output and initialize progress stream
         * @param progress_fn Filename for progress log
         */
        void set_progress_fn (std::string progress_fn) {
            progress_ = std::make_unique<std::ofstream> (progress_fn);
            (*progress_) << "Methods:";
            for (auto str : method_name_)
                (*progress_) << " " << str;
            (*progress_) << std::endl;
        }
        /**
         * @brief Get quiet level
         * @return Quiet level (0=all messages, 1=LNS only, 2=no messages)
         */
        int quiet() const {return quiet_;}
        /**
         * @brief Set quiet level for output control
         * @param quiet Quiet level (0=all messages, 1=LNS only, 2=no messages)
         */
        void set_quiet (int quiet) {
            quiet_ = quiet;
        }
        /**
         * @brief Get current iteration count
         * @return Number of iterations completed
         */
        int iters() const {return iter_;}
        /**
         * @brief Get inner iteration count
         * @return Number of inner iterations
         */
        int inner_iters() const {return inner_iters_;}
        /**
         * @brief Get count of improvements found
         * @return Number of times a better solution was found
         */
        int num_improves() const {return num_improves_;}
        /**
         * @brief Get count of simulated annealing acceptances
         * @return Number of times simulated annealing accepted a worse solution
         */
        int num_simanneal() const {return num_simanneal_;}
        /**
         * @brief Get iteration when best solution was found
         * @return Iteration number of best solution
         */
        int iter_found_best() const {return iter_found_best_;}
        /**
         * @brief Get target flux for an experiment
         * @param exp Experiment index
         * @return Target flux value
         */
        double target_flux(size_t exp) const {return target_flux_[exp];}
        
        /**
         * @brief Write metabolic model to output stream
         * @param out Output stream
         */
        void write_model (std::ostream& out);
        /**
         * @brief Write metabolic model with specific solution to output stream
         * @param out Output stream
         * @param sol Multi-solution to write
         */
        void write_model (std::ostream& out, MultiSolPtr sol);
        /**
         * @brief Write Pareto front to output stream and file
         * @param out Output stream
         * @param save_pareto_fn Filename for Pareto data
         * @param header Header string for output
         */
        void write_pareto (
            std::ostream& out, std::string save_pareto_fn, std::string header
        );
        /**
         * @brief Write solver description to output stream
         * @param out Output stream
         */
        void write_descr (std::ostream& out);

    private:
        /** @brief Sentinel value indicating a reaction is never tabu */
        static constexpr size_t never = 999999;

        /**
         * @enum Method
         * @brief Neighborhood search methods for selecting/adding reactions
         *
         * Defines different heuristics for choosing which reactions to add or remove
         * during the search process. Methods include bias-based, cost-based, flux-based,
         * random selection, and unit cost enforcement.
         */
        enum Method {
            BIAS_SEL, COST_SEL, FLUX_SEL, COST_FLUX_SEL, RAND_SEL, BAD_SEL,
            ADD_COST, ADD_RAND,
            MAKE_UNIT,
            NUM_METHODS
        };
        /** @brief Names of neighborhood search methods for display */
        std::vector<const char*> method_name_ = {
            "BIAS_SEL", "COST_SEL", "FLUX_SEL", "COST_FLUX_SEL", "RAND_SEL", "BAD_SEL",
            "ADD_COST", "ADD_RAND",
            "MAKE_UNIT",
            "NUM_METHODS"
        };
        enum {NUM_SELECT = 5, /**< Number of reaction selection methods */
              NUM_ADD = 2,    /**< Number of reaction addition methods */
              NUM_UNIT = 1    /**< Number of unit cost methods */
        };
        
        /** @brief Shared pointer type for reaction lists */
        using ReactListPtr = std::shared_ptr<ReactList>;

        /**
         * @brief Get number of experiments in scenario
         * @return Number of experiments
         */
        size_t num_exp() const {return scenario_->num_experiments();}

        /**
         * @brief Attempt to improve current solution
         * @return Pointer to improved solution, or nullptr
         */
        SolutionPtr improve ();
        /**
         * @brief Initialize solver state before main solve loop
         */
        void init_solve();
        /**
         * @brief Set all reactions to unit cost
         */
        void make_all_unit_cost();
        /**
         * @brief Unselect reactions that are not gene-independent
         * @param num_unselected Output parameter for count of unselected reactions
         */
        void unselect_non_gene_ind(int& num_unselected);
        /**
         * @brief Solve FBA for all experiments
         * @param first_time Whether this is the first solve
         * @return Pointer to multi-solution
         */
        MultiSolPtr solve_all(bool first_time);
        /**
         * @brief Execute FBA solve and check solution validity
         * @return Pointer to multi-solution
         */
        MultiSolPtr do_solve();
        /**
         * @brief Enable (select) a reaction
         * @param react Pointer to reaction to enable
         */
        void enable_reaction (Reaction* react);
        /**
         * @brief Enable (select) multiple reactions
         * @param react_list List of reactions to enable
         */
        void enable_reactions (ReactList& react_list)
        {
            for (auto react : react_list)
                enable_reaction (react);
        }
        /**
         * @brief Disable (unselect) a reaction
         * @param react Pointer to reaction to disable
         */
        void disable_reaction (Reaction* react);
        /**
         * @brief Disable (unselect) multiple reactions
         * @param react_list List of reactions to disable
         */
        void disable_reactions (ReactList& react_list)
        {
            for (auto react : react_list)
                disable_reaction (react);
        }
        /**
         * @brief Set cost coefficient for a reaction
         * @param react Pointer to reaction
         * @param cost Cost value to set
         */
        void set_react_cost (Reaction* react, double cost);
        
        /**
         * @brief Set reactions to unit cost
         * @param react_list List of reactions to modify
         */
        void set_unit_react_cost (ReactList react_list)
        {
            for (auto react : react_list) {
                set_react_cost (react, unit_cost_);
                non_unit_react_[react->index()] = false;
            }
        }
        /**
         * @brief Reset reactions to their original cost coefficients
         * @param react_list List of reactions to reset
         */
        void reset_react_cost (ReactList react_list)
        {
            for (auto react : react_list) {
                set_react_cost (react, react->obj_coeff());
                non_unit_react_[react->index()] = true;
            }
        }
        /**
         * @brief Reset reaction costs to match target non-unit flags
         * @param match_non_unit Boolean vector indicating target non-unit status
         */
        void reset_react_cost (ReactBools& match_non_unit)
        {
            for (auto& react_ptr : scenario_->reactions()) {
                auto react = react_ptr.get();
                size_t k = react->index();
                if (non_unit_react_[k] != match_non_unit[k]) {
                    if (non_unit_react_[k]) {
                        set_react_cost (react, unit_cost_);
                        non_unit_react_[k] = false;
                    }
                    else {
                        if (!limeDblEqual (react->obj_coeff(), unit_cost_)) {
                            set_react_cost (react, react->obj_coeff());
                            non_unit_react_[k] = true;
                        }
                    }
                }
            }
        }
        
        /**
         * @brief Save solution as best if it improves objectives
         * @param sol Solution attributes to potentially save as best
         */
        void save_best (SolAttributesPtr sol);
        /**
         * @brief Save list of unselected reactions from best solution
         */
        void save_best_unselected ();
        /**
         * @brief Check if a reaction is in a list
         * @param the_list List to search
         * @param react Reaction to find
         * @return True if reaction is in list
         */
        bool list_contains (ReactList& the_list, Reaction* react) {
            return
                std::find(the_list.begin(), the_list.end(), react) !=
                the_list.end();
        }
        /**
         * @brief Select reactions to remove using specified method
         * @param num_to_remove Number of reactions to select
         * @param method Selection method to use
         * @param reacts Output list of selected reactions
         * @return True if selection succeeded
         */
        bool select_reactions_to_remove (
            int num_to_remove, Method method, ReactList& reacts
        );
        /**
         * @brief Select reactions to add using specified method
         * @param num_to_add Number of reactions to select
         * @param method Selection method to use
         * @param reacts Output list of selected reactions
         * @return True if selection succeeded
         */
        bool select_reactions_to_add (
            int num_to_add, Method method, ReactList& reacts
        );
        /**
         * @brief Select cheap reactions to add from a candidate list
         * @param reacts_to_add Candidate reactions to consider
         * @param num_to_add Number of reactions to select
         * @param reacts Output list of selected reactions
         * @return True if selection succeeded
         */
        bool select_cheap_reactions_to_add (
            ReactList& reacts_to_add,
            int num_to_add, ReactList& reacts
        );
        /**
         * @brief Select reactions to convert to unit cost using specified method
         * @param num_to_change Number of reactions to select
         * @param method Selection method to use
         * @param reacts Output list of selected reactions
         * @return True if selection succeeded
         */
        bool select_reactions_to_make_unit (
            int num_to_change, Method method, ReactList& reacts
        );
        /**
         * @brief Select unused reactions to convert to unit cost
         * @param num_to_change Number of reactions to select
         * @param reacts Output list of selected reactions
         * @return True if selection succeeded
         */
        bool select_unused_reactions_to_make_unit (
            int num_to_change, ReactList& reacts
        );
        /**
         * @brief Remove reactions using bias toward target experiment
         * @param target_exp Index of experiment to bias toward
         */
        void remove_reactions_biased (size_t target_exp);
        /**
         * @brief Remove reactions based on cost
         */
        void remove_reactions_cost ();
        /**
         * @brief Remove reactions based on flux values
         */
        void remove_reactions_flux ();
        /**
         * @brief Remove reactions based on combined cost and flux
         */
        void remove_reactions_cost_flux ();
        /**
         * @brief Remove reactions randomly
         */
        void remove_reactions_rand ();
        /**
         * @brief Remove reactions marked as bad
         */
        void remove_reactions_bad ();
        
        /**
         * @brief Check if biomass flux exceeds target
         * @param exp Experiment index
         * @param sol Multi-solution to check
         * @return True if flux is over target
         */
        bool is_over_target (size_t exp, MultiSolPtr& sol)
        {
            return
                sol->biomass_flux(exp) >
                target_flux_[exp] * params_->rank_tol_mult;
        }
        /**
         * @brief Check if biomass flux is on target (not over)
         * @param exp Experiment index
         * @param sol Multi-solution to check
         * @return True if flux is on or below target
         */
        bool is_on_target (size_t exp, MultiSolPtr& sol)
        {
            return !is_over_target (exp, sol);
        }
        /**
         * @brief Check if biomass flux is under target
         * @param exp Experiment index
         * @param sol Multi-solution to check
         * @return True if flux is significantly under target
         */
        bool is_under_target (size_t exp, MultiSolPtr& sol)
        {
            return
                sol->biomass_flux(exp) <
                target_flux_[exp] / params_->rank_tol_mult;
        }
        /**
         * @brief Get number of iterations without improvement
         * @return Iterations since best solution was found
         */
        int iters_no_improve() const {return iter_ - iter_found_best_;}

        /**
         * @brief Check if a reaction is currently tabu
         * @param react Reaction to check
         * @return True if reaction is tabu
         */
        bool is_tabu(const Reaction* react)
        {
            return tabu_list_[react->index()] != never;
        }
        /**
         * @brief Mark reactions as tabu until specified iteration
         * @param the_list List of reactions to make tabu
         * @param iter Iteration when tabu status expires
         */
        void make_tabu (ReactList& the_list, size_t iter)
        {
            for (auto react : the_list) {
                DEBUG (
                    'y', "  React " << *react <<
                    " tabu until iter " << iter
                );
                tabu_list_[react->index()] = iter;
            }
        }
        /**
         * @brief Clear tabu status for a reaction
         * @param react Reaction to clear
         */
        void clear_tabu (const Reaction* react)
        {
            tabu_list_[react->index()] = never;
        }
        /**
         * @brief Check and update tabu list for current iteration
         * @param iter Current iteration number
         */
        void check_tabu_list(size_t iter)
        {
            for (size_t k = 0; k < num_reactions(); k++) {
                if (tabu_list_[k] <= iter) {
                    DEBUG (
                        'y', "  React " << scenario_->reaction(k)->name() <<
                        " no longer tabu"
                    );
                    tabu_list_[k] = never;
                }
            }
        }
        /**
         * @brief Attempt to remove tabu restriction
         * @return True if tabu was removed
         */
        bool remove_tabu ();
        
        /**
         * @brief Position cursor on screen for display
         * @param row Row position
         * @param col Column position
         */
        void screen_pos (int row, int col) {
            if (quiet_ < 2) {
                std::cout << escSeqCmd (lime::POS, row, col);
            }
        }
        /**
         * @brief Clear terminal screen
         * @param whole_screen If true, clear entire screen; otherwise clear from cursor
         */
        void clear_screen (bool whole_screen = true) {
            if (quiet_ < 2) {
                std::cout << escSeqCmd (
                    whole_screen ? lime::CLRSCR_HOME : lime::CLR_REST_SCREEN
                );
            }
        }
        /**
         * @brief Display message on screen without newline
         * @param message Message to display
         */
        void on_screen_nl (std::string message) {
            DEBUG ('A', message);
            if (quiet_ < 2) {
                std::cout <<
                    escSeqColour (lime::CYAN) <<
                    message <<
                    escSeqColour (lime::RESET) <<
                    escSeqCmd (lime::CLR_REST_LINE);
            }
        }
        /**
         * @brief Display message on screen with newline
         * @param message Message to display
         */
        void on_screen (std::string message) {
            DEBUG ('A', message);
            if (quiet_ < 2) {
                std::cout <<
                    escSeqColour (lime::CYAN) <<
                    message <<
                    escSeqColour (lime::RESET) <<
                    escSeqCmd (lime::CLR_REST_LINE) <<
                    std::endl;
            }
        }

        Flavour flavour_; /**< Solution flavor (minimize/maximize biomass) */
        WhichObj which_obj_; /**< Active objective function combination */
        // The elements of the multi-objective optimisation that
        // forms the pareto front - in order
        std::vector<ObjElt> obj_elt_; /**< Elements of multi-objective in order */
#ifdef PLUM_GUROBI
        GrbEnvPtr grbenv_; /**< Gurobi environment pointer (if PLUM_GUROBI defined) */
#else
        void* grbenv_; /**< Placeholder for Gurobi environment (not compiled with Gurobi) */
#endif
        std::vector<LPSolverImpPtr> imp_; /**< LP solver implementations for each experiment */
        std::vector<LPSolverPtr> lp_solver_; /**< LP solver wrappers for each experiment */
        lime::SimAnneal<double> simanneal_; /**< Simulated annealing controller */
        lime::ParetoFront<SolAttributesPtr> pareto_front_; /**< Pareto front of non-dominated solutions */
        
        SolAttributesPtr best_sol_; /**< Best solution found so far */
        SolAttributesPtr incumb_sol_; /**< Current incumbent solution */
        
        std::string progname_; /**< Program name for logging */
        std::string save_best_fn_; /**< Filename for saving best solution */
        std::string save_best_sol_fn_; /**< Filename for saving best solution details */

        std::unique_ptr<std::ofstream> progress_; /**< Output stream for progress logging */
        
        int max_iters_; /**< Maximum number of iterations allowed */
        bool make_all_unit_; /**< Whether to force all reactions to unit cost */
        // Quiet = 0 - lnsmx + gurobi messages
        //         1 - lnsmx messages only
        //         2 - no messages
        int quiet_; /**< Quiet level (0=all messages, 1=LNS only, 2=none) */

        // All solves use a fixed obj mult
        double biomass_obj_mult_; /**< Fixed objective multiplier for biomass */
        // Cost base, set after initial solves
        double base_cost_; /**< Base cost for normalization, set after initial solves */

        int iter_; /**< Current iteration count */
        int inner_iters_; /**< Inner iteration count */
        int num_improves_; /**< Count of improvements found */
        int num_simanneal_; /**< Count of simulated annealing acceptances */
        int iter_found_best_; /**< Iteration when best solution was found */

        std::vector<double> target_flux_; /**< Target flux values for each experiment */

        lime::BiasChoiceT<Reaction*> react_chooser_; /**< Biased random chooser for reactions */
        lime::BiasChoiceT<size_t> exp_chooser_; /**< Biased random chooser for experiments */

        // Adapt the choice of select and add methods
        lime::AdaptChoice adapt_choice_; /**< Adaptive choice controller for search methods */

        // Reactions we don't want to delete. Holds iter when no longer tabu
        std::vector<size_t> tabu_list_; /**< Tabu list holding iteration when each reaction becomes available */

        // The change rank0 target at each interval
        double target_delta_; /**< Change in rank0 target at each interval */
        int interval_len_; /**< Length of interval for target changes */

        // Value used for "unit" costs
        double unit_cost_; /**< Value used for unit costs */
        ReactBools non_unit_react_; /**< Boolean flags indicating non-unit cost reactions */
        // List of reactions we want to get rid of
        ReactList bad_reacts_; /**< List of reactions to preferentially remove */

        // The list of reactions we could not make unit cost,
        // and had to eliminate
        ReactList unit_fails_; /**< Reactions that could not be made unit cost and were eliminated */

        // Kludgy returns from do_solve, only used for logging/debugging
        /**
         * @brief Enumeration of solution failure reasons
         *
         * Used for logging and debugging to track why a solve attempt failed.
         */
        enum {
            NO_FAIL = 0, NO_SOL = 1, RUNAWAYS = 2, NO_GROWTH0 = 3,
            EXCESS_GROWTH0 = 4, GROWTH_MISMATCH = 5
        } why_failed_;
        double rank0_flux_; /**< Rank 0 flux value from last solve (for debugging) */
    };
}
