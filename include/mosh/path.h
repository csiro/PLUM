/**
 * @file path.h
 * @brief Defines graph path structures for metabolic network traversal.
 *
 * This file contains classes representing edges and paths in metabolic networks,
 * used for gap-filling analysis and pathway reconstruction in flux balance analysis.
 */
#pragma once

#include <vector>
#include <memory>


#include "mosh/metabolite.h"
#include "mosh/reaction.h"


/**
 * @brief Namespace for metabolic optimization and simulation heuristics.
 */
namespace mosh
{
    /**
     * @class Edge
     * @brief Represents a directed edge in a metabolic network graph.
     *
     * An Edge connects two metabolites through a specific reaction, forming
     * the basic unit of metabolic pathways. Used in gap-filling algorithms
     * to trace metabolite transformations through the reaction network.
     */
    class Edge {
        
    public:
        /**
         * @brief Constructs an Edge connecting two metabolites via a reaction.
         * @param from Pointer to the source metabolite.
         * @param to Pointer to the destination metabolite.
         * @param via Pointer to the reaction connecting the metabolites.
         */
        Edge (
            Metabolite* from, Metabolite* to, Reaction* via
        ) :
            from_(from),
            to_(to),
            via_(via)
        {
        }
        
    private:
        Metabolite* from_; /**< Source metabolite of the edge. */
        Metabolite* to_; /**< Destination metabolite of the edge. */
        Reaction* via_; /**< Reaction that connects the source to destination metabolite. */"
    };

    /**
     * @typedef EdgePtr
     * @brief Shared pointer type for Edge objects.
     *
     * Used for safe memory management of edges in metabolic pathway graphs.
     */
    using EdgePtr = std::shared_ptr<Edge>;

    /**
     * @class Path
     * @brief Represents a sequence of edges forming a metabolic pathway.
     *
     * A Path is a container of edges that represents a complete route through
     * the metabolic network. Used in gap-filling algorithms to identify and
     * evaluate potential metabolic pathways connecting reactants to products.
     */
    class Path {
        
    public:
        /**
         * @brief Default constructor creating an empty path.
         */
        Path () :
            path_()
        {
        }

        /**
         * @brief Returns an iterator to the beginning of the path.
         * @return Iterator pointing to the first edge in the path.
         */
        std::vector<EdgePtr>::iterator begin() {return path_.begin();}
        /**
         * @brief Returns an iterator to the end of the path.
         * @return Iterator pointing past the last edge in the path.
         */
        std::vector<EdgePtr>::iterator end() {return path_.end();}
        
    private:
        std::vector<EdgePtr> path_; /**< Ordered sequence of edges forming the pathway. */"
        
    };

}
