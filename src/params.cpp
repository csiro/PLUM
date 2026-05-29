
#include "lime/numutil.h"
#include "lime/debug.h"

#include "mosh/params.h"
#include "mosh/reaction.h"

using namespace std;
using namespace lime;
using namespace mosh;

bool
Params::is_gene_indicated (const Reaction* react) const
{
    return
        limeLessEq (react->obj_coeff(), gene_ind_cost) &&
        react->obj_coeff() >= 0.0f;
}

bool
Params::exceeds_max_cost (const Reaction* react) const
{
    return max_react_cost > 0 && react->obj_coeff() > max_react_cost;
}
