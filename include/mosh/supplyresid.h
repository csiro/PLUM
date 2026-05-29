#pragma once

#include <memory>

namespace mosh
{
    class SupplyResid 
    {
    public:
        SupplyResid(double lb, double ub) :
            lb_(lb),
            ub_(ub)
        {
        }

        size_t index() const {return index_;}
        void set_index (size_t index) {index_ = index;}
        bool is_supply() const {return limeIsNegative(sr_lb_);}
        double supply() const {
            assert (is_supply());
            return -sr_lb_;
        }
        bool is_residual() const {return limeIsPositive(sr_ub_);}
        double residual() const {
            assert (is_residual());
            return sr_ub_;
        }
        bool is_balanced() const {return limeIsZero(sr_lb_);}
        double sr_lb() const {return sr_lb_;}
        double sr_ub() const {return sr_ub_;}
        // When set, this value has human-centred meaning:
        // +ve -> this is a supplied metabolite
        // -ve -> this is a residual metabolite
        //   0 -> this is a balanced metabolite
        void set_supply_resid (double lb, double ub) {
            // Negate to save in the mathematically-centred meaning
            sr_lb_ = limeMin (-lb, -ub);
            sr_ub_ = limeMax (-lb, -ub);
            DEBUG (
                'm', "  Set lb/ub to " << sr_lb_ << " " << sr_ub_ <<
                " " << sr_ub() 
            );
        }
    
    private:
        std::string name_;
        size_t index_;
        // Supply/residual upper and lower bound
        // -ve -> this is a supplied metabolite
        // +ve -> this is a residual metabolite
        //   0 -> this is a balanced metabolite
        double sr_lb_;  
        double sr_ub_;
        double dummy_;
    };
    
    using SupplyResidPtr = std::shared_ptr<SupplyResid>;
}
