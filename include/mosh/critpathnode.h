#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <memory>

#include "lime/spnode.h"
#include "lime/numutil.h"

#include "mosh/reaction.h"

namespace mosh
{
    class CritPathNode : public lime::SPNode
    {
    public:
        CritPathNode (const Reaction* react, int cost) :
            SPNode(cost),
            react_(react)
        {
        }

        const Reaction* react() const {return react_;}
        
        bool isLessThan(const SPNode* other_sp) const override
        {
            const CritPathNode* other = (const CritPathNode*) other_sp;
            return react_->index() < other->react_->index();
        }
    
    private:
        const Reaction* react_;
    };
    
    using ReactionPtr = std::shared_ptr<Reaction>;
}
