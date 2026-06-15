/**
 * @file dualvals.h
 * @brief Storage container for dual values from linear programming solutions
 *
 * This file defines the DualVals class which stores lower bound and upper bound
 * dual values (shadow prices) for metabolite constraints obtained from flux balance
 * analysis solutions. Dual values indicate the marginal change in the objective
 * function with respect to changes in metabolite concentration bounds.
 */
#pragma once

/**
 * @brief Main namespace for the MOSH metabolic optimization system
 */
namespace mosh
{
    /**
     * @class DualVals
     * @brief Container for dual values (shadow prices) from LP solutions
     *
     * Stores the lower bound and upper bound dual values for each metabolite
     * constraint in a flux balance analysis problem. Dual values represent the
     * sensitivity of the objective function to changes in metabolite bounds and
     * are used in gap-filling algorithms to identify important metabolic constraints.
     */
    class DualVals
    {
    public:
        /**
         * @brief Constructor that initializes dual value storage for a given number of metabolites
         * @param num_metabolites Number of metabolites in the metabolic network
         */
        DualVals (size_t num_metabolites) :
            lb_dual_(num_metabolites),
            ub_dual_(num_metabolites)
        {
        }
            
        /**
         * @brief Get the lower bound dual value for a metabolite by index
         * @param k Index of the metabolite
         * @return Lower bound dual value (shadow price) for the metabolite constraint
         */
        double lb_dual(size_t k) const {return lb_dual_[k];}
        /**
         * @brief Get the upper bound dual value for a metabolite by index
         * @param k Index of the metabolite
         * @return Upper bound dual value (shadow price) for the metabolite constraint
         */
        double ub_dual(size_t k) const {return ub_dual_[k];}
        /**
         * @brief Get the lower bound dual value for a metabolite by pointer
         * @param met Pointer to the Metabolite object
         * @return Lower bound dual value (shadow price) for the metabolite constraint
         */
        double lb_dual(const Metabolite* met) const {
            return lb_dual (met->index());
        }
        /**
         * @brief Get the upper bound dual value for a metabolite by pointer
         * @param met Pointer to the Metabolite object
         * @return Upper bound dual value (shadow price) for the metabolite constraint
         */
        double ub_dual(const Metabolite* met) const {
            return ub_dual (met->index());
        }
        /**
         * @brief Set both lower and upper bound dual values for a metabolite
         * @param k Index of the metabolite
         * @param lb_dual Lower bound dual value to set
         * @param ub_dual Upper bound dual value to set
         */
        void set_duals (size_t k, double lb_dual, double ub_dual) {
            lb_dual_[k] = lb_dual;
            ub_dual_[k] = ub_dual;
        }
        
    private:
        /**
         * @brief Default constructor (private, not intended for external use)
         */
        DualVals () {}
        
        std::vector<double> lb_dual_; /**< Lower bound dual values for each metabolite */
        std::vector<double> ub_dual_; /**< Upper bound dual values for each metabolite */
    };
    /**
     * @typedef DualValsPtr
     * @brief Shared pointer type for DualVals objects
     */
    using DualValsPtr = std::shared_ptr<DualVals>;
}
