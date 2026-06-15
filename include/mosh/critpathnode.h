/**
 * @file critpathnode.h
 * @brief Defines the CritPathNode class for representing reactions in critical path computations.
 *
 * This file contains the CritPathNode class which extends lime::SPNode to represent
 * reactions in shortest path algorithms used for metabolic gap-filling and flux balance
 * analysis. Each node represents a reaction with an associated cost in the critical path.
 */

#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <memory>

#include "lime/spnode.h"
#include "lime/numutil.h"

#include "mosh/reaction.h"

/**
 * @namespace mosh
 * @brief Namespace for metabolic optimization and shortest path heuristics.
 */
namespace mosh
{
    /**
     * @class CritPathNode
     * @brief Represents a reaction node in critical path computations for metabolic networks.
     *
     * CritPathNode extends lime::SPNode to provide a specialized node type for shortest
     * path algorithms applied to metabolic networks. Each node wraps a Reaction object
     * and maintains an associated cost for use in gap-filling and optimization algorithms.
     * The nodes can be compared based on their underlying reaction indices to maintain
     * consistent ordering in priority queues and other data structures.
     */
    class CritPathNode : public lime::SPNode
    {
    public:
        /**
         * @brief Constructs a CritPathNode with a reaction and associated cost.
         * @param react Pointer to the Reaction object this node represents (must not be null).
         * @param cost The cost associated with including this reaction in the path.
         */
        CritPathNode (const Reaction* react, int cost) :
            SPNode(cost),
            react_(react)
        {
        }

        /**
         * @brief Returns the reaction associated with this node.
         * @return Const pointer to the Reaction object.
         */
        const Reaction* react() const {return react_;}

        /**
         * @brief Compares this node with another based on reaction indices.
         * @param other_sp Pointer to another SPNode (expected to be a CritPathNode).
         * @return true if this node's reaction index is less than the other's, false otherwise.
         *
         * This comparison is used to maintain consistent ordering of nodes in data structures
         * such as priority queues during shortest path computations.
         */
        bool isLessThan(const SPNode* other_sp) const override
        {
            const CritPathNode* other = (const CritPathNode*) other_sp;
            return react_->index() < other->react_->index();
        }
    
    private:
        const Reaction* react_; /**< Pointer to the reaction represented by this node. */"
    };
    
    /**
     * @typedef ReactionPtr
     * @brief Shared pointer type for managing Reaction object lifetimes.
     */
    using ReactionPtr = std::shared_ptr<Reaction>;
}
