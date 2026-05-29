#pragma once

namespace mosh
{
    class DualVals
    {
    public:
        DualVals (size_t num_metabolites) :
            lb_dual_(num_metabolites),
            ub_dual_(num_metabolites)
        {
        }
            
        double lb_dual(size_t k) const {return lb_dual_[k];}
        double ub_dual(size_t k) const {return ub_dual_[k];}
        double lb_dual(const Metabolite* met) const {
            return lb_dual (met->index());
        }
        double ub_dual(const Metabolite* met) const {
            return ub_dual (met->index());
        }
        void set_duals (size_t k, double lb_dual, double ub_dual) {
            lb_dual_[k] = lb_dual;
            ub_dual_[k] = ub_dual;
        }
        
    private:
        DualVals () {}
        
        std::vector<double> lb_dual_;
        std::vector<double> ub_dual_;
    };
    using DualValsPtr = std::shared_ptr<DualVals>;
}
