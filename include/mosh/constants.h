#pragma once

namespace mosh
{
    // Cost of dummy reactions
    constexpr double DUMMY_FLUX_UB = 1000.0f;
    
    constexpr int IGNORE = 9999;
    constexpr int DONT_CARE = 9999;

    bool isDontCare (double val);
    bool isIgnoreVal (double val);

    enum StatusEnum {
        UNSOLVED, CTS_OPTIMAL, INT_OPTIMAL, SUB_OPTIMAL, INFEASIBLE_, RANDOM
    };
    enum Flavour {GUROBI, LP_SOLVE, HIGHS};
}
