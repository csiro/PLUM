
#include "lime/debug.h"
#include "lime/linereader.h"
#include "lime/limetok.h"
#include "lime/numutil.h"

#include "mosh/dummysolver.h"

using namespace std;
using namespace lime;
using namespace mosh;


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

