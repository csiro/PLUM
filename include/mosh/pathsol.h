#pragma once

#include <vector>
#include <memory>
#include <set>

#include "mosh/scenario.h"
#include "mosh/reaction.h"
#include "mosh/path.h"

#include "lime/dijkstra.h"

namespace mosh
{
    class PathSol 
    {
    public:
        PathSol(const Scenario* scenario) :
            scenario_(scenario)
        {
        }
        PathSol(PathSol* other) :
            scenario_(other->scenario_)
        {
        }

        const Scenario* scenario() const {return scenario_;}
        
        bool calc_path (Metabolite* from, Metabolite* to, Path& path);

    private:
        lime::Dijkstra<int> dijkstra_;
        
        const Scenario* scenario_;
        Path path_;
    };

    using PathSolPtr = std::shared_ptr<PathSol>;
}
