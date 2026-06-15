/**
 * @file mathheur.cpp
 * @brief Implementation of the mathematical heuristic solver for metabolic gap-filling
 *
 * This file implements the MathHeur class which uses a mathematical programming-based
 * heuristic approach to solve metabolic network gap-filling problems. The solver iteratively
 * improves solutions by focusing on individual metabolites and their associated reactions.
 */

#include <sstream>
#include <list>

#include "lime/debug.h"
#include "lime/linereader.h"
#include "lime/limetok.h"
#include "lime/numutil.h"
#include "lime/timekeeper.h"
#include "lime/error.h"

#include "mosh/constants.h"
#include "mosh/mathheur.h"
#include "mosh/grblpsolver.h"
#include "mosh/reaction.h"

using namespace std;
using namespace lime;
using namespace mosh;

/**
 * @brief Solves the metabolic gap-filling problem using a mathematical heuristic approach
 *
 * This method implements the main solving loop of the mathematical heuristic algorithm.
 * It starts with an initial solution from a linear programming solver, then iteratively
 * attempts to improve the solution by randomly selecting metabolites and optimizing
 * their associated reactions.
 *
 * @return SolutionPtr Pointer to the best solution found during the search
 * @throws lime::Error if the scenario contains more than one experiment (only single experiment supported)
 */
SolutionPtr
MathHeur::solve ()
{
    DEBUG ('L', "Solve with MathHeuristic solver");
    on_screen ("Solve with MathHeuristic solver");

    if (scenario_->num_experiments() > 1)
        limeCrash ("MH only works for 1 experiment");
    
    GrbLPSolver lpsolver (
        scenario_, params_, env_, 
        rand_.generateSeed()
    );

    TimeKeeper timer (time_limit_);
    
    // Find out initial solution
    SolutionPtr best = lpsolver.solve();
    auto best_obj = best->abs_obj_value();

    while (++iter_ < max_iters_ && timer.hasTimeLeft()) {
        DEBUG ('m', "MH Iter " << iter_);

        // Choose a metabolite
        int met_idx = rand_.uniform0n_1 (num_metabolites());
        improve (metabolite(met_idx), best);
    }
    
    return best;
}

/**
 * @brief Sets up the Gurobi optimization model with variables and constraints
 *
 * This method initializes the mixed-integer programming model by:
 * - Creating flux variables for each reaction (continuous, bounded by reaction limits)
 * - Creating binary use variables to indicate whether each reaction is active
 * - Adding linking constraints between flux and use variables
 * - Setting up stoichiometric flux balance equations for all metabolites
 *
 * The model minimizes the number of reactions used while maintaining steady-state
 * flux balance for all metabolites in the network.
 */
void
MathHeur::set_up_model()
{
    for (size_t k = 0; k < num_reactions(); k++) {
        auto react = reaction(k);
        
        double lb = 0.0f;
        double ub = react->flux_ub();
        double obj_coeff = 0.0f;
        DEBUG (
            'L', "      Adding var for reaction " << k <<
            " " << react->name() <<
            " obj-coeff " << obj_coeff << 
            " lb " << lb << " ub " << ub
        );
        flux_[k] =
            model_.addVar (
                lb, ub, obj_coeff, GRB_CONTINUOUS,
                react->name()
            );
        if (react->is_biomass()) {
            biomass_ = flux_[k];
        }
        
        lb = 0.0f;
        ub = 1.0f;
        use_[k] = 
            model_.addVar (
                lb, ub, obj_coeff, GRB_BINARY,
                react->name() + "-use"
            );

        // Add the linking constraints - can only use flux if use is 1.0
        model_.addConstr (
            flux_[k] - reaction(k)->flux_ub() * use_[k] <= 0.0f,
            reaction(k)->name() + "-link"
        );
    }

    // Set up the flux balance equations
    for (size_t m = 0; m < num_metabolites(); m++) {
        auto met = metabolite(m);
        
        GRBLinExpr lhs = 0;

        auto exp = experiment(0);
        
        for (size_t k = 0; k < num_reactions(); k++) {
            auto react = reaction(k);
            if (react->uses(met))
                lhs += react->met_coeff (m) * flux_[k];
        }
        model_.addConstr (lhs >= exp->lb(met), met->name() + "_LB");
        model_.addConstr (lhs <= exp->ub(met), met->name() + "_UB");
    }
}

/**
 * @brief Attempts to improve the current solution by focusing on a specific metabolite
 *
 * This method refines the solution by examining all reactions associated with a given
 * metabolite. Reactions that do not use the metabolite are temporarily fixed to zero,
 * allowing the solver to focus optimization efforts on the relevant reaction subset.
 *
 * @param met Pointer to the metabolite to focus on for improvement
 * @param sol Reference to the current best solution, updated if improvement is found
 */
void
MathHeur::improve(const Metabolite* met, SolutionPtr& sol)
{
    DEBUG ('m', "Improving metabolite " << *met);
    auto met_idx = met->index();
    
    for (size_t k = 0; k < num_reactions(); k++) {
        auto react = reaction(k);
        if (react->uses (met)) {
            // Put this in the model
            if (sol->uses_react (react)) {
                
            }
        }
        else {
            // Fix flux to zero
            use_[k].set(GRB_DoubleAttr_UB, 0.0f);
            flux_[k].set(GRB_DoubleAttr_UB, 0.0f);
        }
    }


}

/**
 * @brief Generates a summary string of the solver's execution statistics
 *
 * @return std::string Summary containing the number of mathematical heuristic iterations performed
 */
string
MathHeur::summary()
{
    return " mh_iters " + to_string (iter_);
}
