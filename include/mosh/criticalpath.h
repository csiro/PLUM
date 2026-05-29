#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <memory>

#include "lime/dijkstra.h"
#include "lime/numutil.h"
#include "lime/dijkstra.h"

#include "mosh/scenario.h"
#include "mosh/solution.h"
#include "mosh/reaction.h"

namespace mosh
{
    class CriticalPath 
    {
    public:
        CriticalPath (
            SolutionPtr sol,
            const Metabolite* carbon_source
        ) :
            sol_(sol),
            scenario_(sol->scenario()),
            carbon_source_(carbon_source),
            graph_(scenario_->num_metabolites())
        {
            init_graph();
        }
            
        void find_paths(std::ostream& out);

    private:
        void init_graph();
        
        SolutionPtr sol_;
        const Scenario* scenario_;
        const Metabolite* carbon_source_;
        lime::Dijkstra<double> graph_;
    };
}
