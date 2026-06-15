/**
 * @file params.cpp
 * @brief Implementation of algorithm parameter utilities for metabolic gap-filling
 *
 * This file contains member function implementations for the Params class,
 * which provides parameter management and validation methods used in flux
 * balance analysis and metabolic network gap-filling algorithms.
 */

#include "lime/numutil.h"
#include "lime/debug.h"

#include "mosh/params.h"
#include "mosh/reaction.h"

using namespace std;
using namespace lime;
using namespace mosh;

/**
 * @brief Determines if a reaction is gene-indicated based on its objective coefficient
 *
 * A reaction is considered gene-indicated if its objective coefficient (cost) is
 * less than or equal to the gene indication cost threshold and is non-negative.
 * Gene-indicated reactions represent pathways that have genomic evidence supporting
 * their presence in the metabolic network.
 *
 * @param react Pointer to the reaction to evaluate
 * @return true if the reaction's cost indicates genomic evidence, false otherwise
 */
bool
Params::is_gene_indicated (const Reaction* react) const
{
    return
        limeLessEq (react->obj_coeff(), gene_ind_cost) &&
        react->obj_coeff() >= 0.0f;
}

/**
 * @brief Checks if a reaction's cost exceeds the maximum allowed cost threshold
 *
 * This method evaluates whether a reaction's objective coefficient surpasses the
 * configured maximum reaction cost. Reactions exceeding this threshold may be
 * excluded from consideration during metabolic gap-filling to limit the search
 * space to biologically plausible solutions. A max_react_cost value less than or
 * equal to zero disables this check.
 *
 * @param react Pointer to the reaction to evaluate
 * @return true if max_react_cost is positive and the reaction cost exceeds it, false otherwise
 */
bool
Params::exceeds_max_cost (const Reaction* react) const
{
    return max_react_cost > 0 && react->obj_coeff() > max_react_cost;
}
