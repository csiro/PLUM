#pragma once

#include <vector>
#include <memory>


#include "mosh/metabolite.h"
#include "mosh/reaction.h"


namespace mosh
{
    class Edge {
        
    public:
        Edge (
            Metabolite* from, Metabolite* to, Reaction* via
        ) :
            from_(from),
            to_(to),
            via_(via)
        {
        }
        
    private:
        Metabolite* from_;
        Metabolite* to_;
        Reaction* via_;
    };

    using EdgePtr = std::shared_ptr<Edge>;

    class Path {
        
    public:
        Path () :
            path_()
        {
        }

        std::vector<EdgePtr>::iterator begin() {return path_.begin();}
        std::vector<EdgePtr>::iterator end() {return path_.end();}
        
    private:
        std::vector<EdgePtr> path_;
        
    };

}
