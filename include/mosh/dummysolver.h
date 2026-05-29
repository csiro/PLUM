#pragma once

#include <vector>

#include "lime/rand.h"

#include "mosh/gapsolver.h"

namespace mosh
{
    class DummySolver : public GapSolver 
    {
    public:
        DummySolver(Scenario* scenario, Params* params, int seed) :
            GapSolver(scenario, params, seed)
        {
        }

        SolutionPtr solve () override;
        
    private:
        lime::Rand rand_;
    };
}
