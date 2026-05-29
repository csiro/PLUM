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

namespace mosh
{
    class LnsMxSolver : public GapSolver
    {
    public:

        using ReactList = std::list<Reaction*>;
        using ReactBools = std::vector<bool>;
        
        enum ObjElt {COST, TAU, ERROR};
        std::vector<const char*> objelt_name_ = {
            "COST", "TAU", "ERROR"
        };
        
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
        std::vector<const char*> obj_name_ = {
            "ERROR/COST", "TAU+ERROR/COST", "TAU+COST/ERROR", "COST+TAU/ERROR",
            "ERROR+COST/TAU", "TAU+ERROR+COST/COST", "TAU/COST", "COST/TAU"
        };
        
        class SolAttributes
        {
        public:
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

            MultiSolPtr& sol() {return sol_;}
            SolutionPtr sol(size_t k) {return sol_->sol(k);}
            double obj() const {return obj_;}
            double obj2() const {return obj_;}
            std::vector<double>& multi_obj() {return multi_obj_;}
            int growth_mismatch() const {return growth_mismatch_;}
            int pc_true_positive() const {return pc_true_pos_;}
            int pc_true_negative() const {return pc_true_neg_;}
            double error() const {return error_;}
            double tau() const {return tau_;}
            double cost() const {return cost_;}
            double cost_ratio() const {return cost_ratio_;}
            double max_cost() const {return max_cost_;}
            int num_runaways() const {return num_runaways_;}
            double rank0_target() const {return rank0_target_;}
            ReactList& unselected() {return unselected_;}
            size_t num_unselected() {return unselected_.size();}
            ReactBools& non_unit_react() {return non_unit_react_;}
            size_t num_non_unit() const {return num_non_unit_;}
            size_t non_unit_used() const {return non_unit_used_;}

            double biomass_flux () const
            {
                return sol_ == nullptr ? 0.0f : sol_->biomass_flux();
            }

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
            MultiSolPtr sol_;
            double obj_;
            double obj2_;
            std::vector<double> multi_obj_;
            int growth_mismatch_;
            int pc_true_pos_;
            int pc_true_neg_;
            double error_;
            double tau_;
            double cost_;
            double cost_ratio_;
            double max_cost_;
            int num_runaways_;
            double rank0_target_;
            ReactList unselected_;
            ReactBools non_unit_react_;
            size_t num_non_unit_;
            size_t non_unit_used_;
        };
        using SolAttributesPtr = std::shared_ptr<SolAttributes>;        
        
        LnsMxSolver(
            Scenario* scenario, Params* params,
            Flavour flavour, WhichObj which_obj,
            int seed, int max_iters, bool make_all_unit, std::string progname
        );

        virtual SolutionPtr solve() override;
        SolutionPtr solve_by_add();

        WhichObj which_obj() const {return which_obj_;}
        ObjElt obj_elt (size_t k) const {return obj_elt_[k];}
        
        MultiSolPtr best_sol() const
        {
            return best_sol_ == nullptr ? nullptr : best_sol_->sol();
        }
        SolAttributesPtr best_attr() {return best_sol_;}

        double calc_kendall_tau (MultiSolPtr sol) const;
        double mean_squared_error (MultiSolPtr sol) const;
        int num_used () const;
        size_t num_non_unit() const
        {
            return best_sol_ == nullptr ? 0 : best_sol_->num_non_unit();
        }
        size_t non_unit_used() const
        {
            return best_sol_ == nullptr ? 0 : best_sol_->non_unit_used();
        }
        // Percent of true positive growth experiments
        int pc_true_positive() const
        {
            return best_sol_ == nullptr ? 0 : best_sol_->pc_true_positive();
        }
        // Percent of true negative growth experiments
        int pc_true_negative() const
        {
            return best_sol_ == nullptr ? 0 : best_sol_->pc_true_negative();
        }
        
        void write_best_unselected (std::ostream& out);
        void write_best_sol (std::ostream& out);
        
        std::string save_best_fn() const {return save_best_fn_;}
        void set_save_best_fn (std::string save_best_fn) {
            save_best_fn_ = save_best_fn;
        }
        std::string save_best_sol_fn() const {return save_best_sol_fn_;}
        void set_save_best_sol_fn (std::string save_best_sol_fn) {
            save_best_sol_fn_ = save_best_sol_fn;
        }
        void set_progress_fn (std::string progress_fn) {
            progress_ = std::make_unique<std::ofstream> (progress_fn);
            (*progress_) << "Methods:";
            for (auto str : method_name_) 
                (*progress_) << " " << str;
            (*progress_) << std::endl;
        }
        int quiet() const {return quiet_;}
        void set_quiet (int quiet) {
            quiet_ = quiet;
        }
        int iters() const {return iter_;}
        int inner_iters() const {return inner_iters_;}
        int num_improves() const {return num_improves_;}
        int num_simanneal() const {return num_simanneal_;}
        int iter_found_best() const {return iter_found_best_;}
        double target_flux(size_t exp) const {return target_flux_[exp];}
        
        void write_model (std::ostream& out);
        void write_model (std::ostream& out, MultiSolPtr sol);
        void write_pareto (
            std::ostream& out, std::string save_pareto_fn, std::string header
        );
        void write_descr (std::ostream& out);

    private:
        static constexpr size_t never = 999999;

        enum Method {
            BIAS_SEL, COST_SEL, FLUX_SEL, COST_FLUX_SEL, RAND_SEL, BAD_SEL,
            ADD_COST, ADD_RAND, 
            MAKE_UNIT,
            NUM_METHODS
        };
        std::vector<const char*> method_name_ = {
            "BIAS_SEL", "COST_SEL", "FLUX_SEL", "COST_FLUX_SEL", "RAND_SEL", "BAD_SEL", 
            "ADD_COST", "ADD_RAND",
            "MAKE_UNIT",
            "NUM_METHODS"
        };
        enum {NUM_SELECT = 5, NUM_ADD = 2, NUM_UNIT = 1};
        
        using ReactListPtr = std::shared_ptr<ReactList>;

        size_t num_exp() const {return scenario_->num_experiments();}

        SolutionPtr improve ();
        void init_solve();
        void make_all_unit_cost();
        void unselect_non_gene_ind(int& num_unselected);
        MultiSolPtr solve_all(bool first_time);
        MultiSolPtr do_solve();
        void enable_reaction (Reaction* react);
        void enable_reactions (ReactList& react_list)
        {
            for (auto react : react_list)
                enable_reaction (react);
        }
        void disable_reaction (Reaction* react);
        void disable_reactions (ReactList& react_list)
        {
            for (auto react : react_list)
                disable_reaction (react);
        }
        void set_react_cost (Reaction* react, double cost);
        
        void set_unit_react_cost (ReactList react_list)
        {
            for (auto react : react_list) {
                set_react_cost (react, unit_cost_);
                non_unit_react_[react->index()] = false;
            }
        }
        void reset_react_cost (ReactList react_list)
        {
            for (auto react : react_list) {
                set_react_cost (react, react->obj_coeff());
                non_unit_react_[react->index()] = true;
            }
        }
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
        
        void save_best (SolAttributesPtr sol);
        void save_best_unselected ();
        bool list_contains (ReactList& the_list, Reaction* react) {
            return
                std::find(the_list.begin(), the_list.end(), react) !=
                the_list.end();
        }
        bool select_reactions_to_remove (
            int num_to_remove, Method method, ReactList& reacts
        );
        bool select_reactions_to_add (
            int num_to_add, Method method, ReactList& reacts
        );
        bool select_cheap_reactions_to_add (
            ReactList& reacts_to_add,
            int num_to_add, ReactList& reacts
        );
        bool select_reactions_to_make_unit (
            int num_to_change, Method method, ReactList& reacts
        );
        bool select_unused_reactions_to_make_unit (
            int num_to_change, ReactList& reacts
        );
        void remove_reactions_biased (size_t target_exp);
        void remove_reactions_cost ();
        void remove_reactions_flux ();
        void remove_reactions_cost_flux ();
        void remove_reactions_rand ();
        void remove_reactions_bad ();
        
        bool is_over_target (size_t exp, MultiSolPtr& sol)
        {
            return
                sol->biomass_flux(exp) >
                target_flux_[exp] * params_->rank_tol_mult;
        }
        bool is_on_target (size_t exp, MultiSolPtr& sol)
        {
            return !is_over_target (exp, sol);
        }
        bool is_under_target (size_t exp, MultiSolPtr& sol)
        {
            return
                sol->biomass_flux(exp) <
                target_flux_[exp] / params_->rank_tol_mult;
        }
        int iters_no_improve() const {return iter_ - iter_found_best_;}

        bool is_tabu(const Reaction* react)
        {
            return tabu_list_[react->index()] != never;
        }
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
        void clear_tabu (const Reaction* react)
        {
            tabu_list_[react->index()] = never;
        }
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
        bool remove_tabu ();
        
        void screen_pos (int row, int col) {
            if (quiet_ < 2) {
                std::cout << escSeqCmd (lime::POS, row, col);
            }
        }
        void clear_screen (bool whole_screen = true) {
            if (quiet_ < 2) {
                std::cout << escSeqCmd (
                    whole_screen ? lime::CLRSCR_HOME : lime::CLR_REST_SCREEN
                );
            }
        }
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

        Flavour flavour_;
        WhichObj which_obj_;
        // The elements of the multi-objective optimisation that
        // forms the pareto front - in order 
        std::vector<ObjElt> obj_elt_;
#ifdef PLUM_GUROBI
        GrbEnvPtr grbenv_;
#else
        void* grbenv_;
#endif
        std::vector<LPSolverImpPtr> imp_;
        std::vector<LPSolverPtr> lp_solver_;
        lime::SimAnneal<double> simanneal_;
        lime::ParetoFront<SolAttributesPtr> pareto_front_;
        
        SolAttributesPtr best_sol_;
        SolAttributesPtr incumb_sol_;
        
        std::string progname_;
        std::string save_best_fn_;
        std::string save_best_sol_fn_;

        std::unique_ptr<std::ofstream> progress_;
        
        int max_iters_;
        bool make_all_unit_;
        // Quiet = 0 - lnsmx + gurobi messages
        //         1 - lnsmx messages only
        //         2 - no messages 
        int quiet_;

        // All solves use a fixed obj mult 
        double biomass_obj_mult_;
        // Cost base, set after initial solves
        double base_cost_;

        int iter_;
        int inner_iters_;
        int num_improves_;
        int num_simanneal_;
        int iter_found_best_;

        std::vector<double> target_flux_;

        lime::BiasChoiceT<Reaction*> react_chooser_;
        lime::BiasChoiceT<size_t> exp_chooser_;

        // Adapt the choice of select and add methods
        lime::AdaptChoice adapt_choice_;

        // Reactions we don't want to delete. Holds iter when no longer tabu
        std::vector<size_t> tabu_list_;

        // The change rank0 target at each interval
        double target_delta_;
        int interval_len_;

        // Value used for "unit" costs
        double unit_cost_;
        ReactBools non_unit_react_;
        // List of reactions we want to get rid of
        ReactList bad_reacts_;

        // The list of reactions we could not make unit cost,
        // and had to eliminate
        ReactList unit_fails_;

        // Kludgy returns from do_solve, only used for logging/debugging
        enum {
            NO_FAIL = 0, NO_SOL = 1, RUNAWAYS = 2, NO_GROWTH0 = 3,
            EXCESS_GROWTH0 = 4, GROWTH_MISMATCH = 5
        } why_failed_;
        double rank0_flux_;
    };
}
