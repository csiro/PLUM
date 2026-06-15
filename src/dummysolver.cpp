/** @file dummysolver.cpp
 *  @brief Implementation of DummySolver for metabolic network gap-filling with random flux assignment.
 *
 *  This file implements the DummySolver class, which provides a simple randomized solver
 *  for testing purposes. Instead of performing actual flux balance analysis or optimization,
 *  it assigns random flux values to selected reactions in the metabolic network.
 */

#include "lime/debug.h"
#include "lime/linereader.h"
#include "lime/limetok.h"
#include "lime/numutil.h"

#include "mosh/dummysolver.h"

using namespace std;
using namespace lime;
using namespace mosh;


/**
 * @brief Solves the metabolic gap-filling problem by assigning random flux values.
 *
 * This is a dummy implementation that does not perform actual optimization.
 * For each selected reaction in the scenario, it assigns a random flux value
 * between 0 and 2 using a uniform distribution. This solver is primarily used
 * for testing and debugging purposes rather than producing biologically meaningful
 * flux distributions.
 *
 * @return SolutionPtr A shared pointer to a Solution object containing random flux
 *         assignments for all selected reactions. The solver status is set to RANDOM.
 *
 * @note This method does not guarantee mass balance or thermodynamic feasibility.
 * @see GapSolver::solve()
 */
SolutionPtr
DummySolver::solve ()
{
    SolutionPtr sol = make_shared<Solution> (scenario_, params_);

    for (size_t k = 0; k < num_reactions(); k++) {
        if (reaction(k)->is_selected()) {
            // Random val between 0 and 2
            sol->set_flux (k, 2.0f * rand_.uniform01());
        }
    }
    status_ = RANDOM;
    return sol;
}

