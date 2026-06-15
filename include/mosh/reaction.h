/**
 * @file reaction.h
 * @brief Defines the Reaction class for metabolic network analysis
 *
 * This file contains the Reaction class which represents a biochemical reaction
 * in a metabolic network. Reactions are used in flux balance analysis (FBA) and
 * metabolic gap-filling algorithms to model the transformation of metabolites
 * with associated stoichiometric coefficients and flux constraints.
 */

#pragma once

#include <iostream>
#include <assert.h>
#include <vector>
#include <map>
#include <memory>

#include "lime/displayable.h"
#include "lime/numutil.h"
#include "lime/strutil.h"
#include "lime/debug.h"

#include "mosh/metabolite.h"

/**
 * @brief Namespace for metabolic optimization and systems handling
 */
namespace mosh
{
    /**
     * @class Reaction
     * @brief Represents a biochemical reaction in a metabolic network
     *
     * A Reaction models the transformation of metabolites (substrates to products)
     * with stoichiometric coefficients. It supports flux balance analysis by tracking
     * flux bounds, objective coefficients, and reaction properties such as reversibility,
     * exchange reactions (EX), and demand reactions (DM). Reactions can be active/inactive
     * and selected/deselected for linear programming optimization.
     */
    class Reaction : public lime::Displayable
    {
    public:
        /**
         * @brief Constructs a new Reaction object
         *
         * @param name Short name identifier for the reaction
         * @param obj_coeff Objective function coefficient (negative for biomass reactions)
         * @param flux_ub Upper bound for reaction flux
         * @param is_active Whether the reaction is currently active in the network
         * @param is_selected Whether the reaction is selected for LP optimization
         * @param is_dummy Whether this is a dummy reaction (default: false)
         */
        Reaction (
            std::string name, double obj_coeff, double flux_ub,
            bool is_active, bool is_selected,
            bool is_dummy = false
        ) :
            Displayable(),
            name_(name),
            full_name_(name),
            formula_from_(""),
            formula_to_(""),
            obj_coeff_(obj_coeff),
            is_biomass_(limeIsNegative(obj_coeff)),
            is_export_(!is_dummy),
            mets_(),
            in_mets_(),
            out_mets_(),
            met_coeff_(),
            flux_ub_(flux_ub),
            has_known_flux_(false),
            known_flux_(0.0f),
            known_flux_error_(0.0f),
            is_dummy_(is_dummy),
            is_active_(is_active),
            is_selected_(is_selected),
            reduced_cost_(0.0f)
        {
        }

        /**
         * @brief Gets the reaction name
         * @return The short name identifier of the reaction
         */
        std::string name() const override {return name_;}
        /**
         * @brief Sets the reaction name
         * @param name New short name for the reaction
         */
        void set_name(std::string name) {name_ = name;}
        /**
         * @brief Gets the full descriptive name of the reaction
         * @return The full name string
         */
        std::string full_name() const {return full_name_;}
        /**
         * @brief Sets the full descriptive name
         * @param full_name New full name for the reaction
         */
        void set_full_name(std::string full_name) {full_name_ = full_name;}
        /**
         * @brief Gets the reaction formula as a string
         * @return Formula in the format "substrates -> products"
         */
        std::string formula() const
        {
            return formula_from_ + "->" + formula_to_;
        }
        /**
         * @brief Gets a human-readable version of the reaction formula
         * @return Formula with special characters decoded for display
         */
        std::string nice_formula() const
        {
            return nicify (formula_from_) + "->" + nicify (formula_to_);
        }
        /**
         * @brief Converts encoded characters in a name to readable format
         * @param name Name string with encoded special characters
         * @return String with special characters decoded (e.g., __40__ to '(')
         */
        std::string nicify (std::string name) const
        {
            std::string nice_name = name;
            lime::replaceAll (nice_name, "__40__", "(");
            lime::replaceAll (nice_name, "__41__", ")");
            lime::replaceAll (nice_name, "__91__", "[");
            lime::replaceAll (nice_name, "__93__", "]");
            return nice_name;
        }
        /**
         * @brief Gets the reaction index in the network
         * @return Zero-based index of the reaction
         */
        size_t index() const {return index_;}
        /**
         * @brief Sets the reaction index
         * @param index New index value for the reaction
         */
        void set_index (size_t index) {index_ = index;}
        /**
         * @brief Sets the objective function coefficient
         * @param obj_coeff New objective coefficient value
         */
        void set_obj_coeff (double obj_coeff) {obj_coeff_ = obj_coeff;}
        /**
         * @brief Gets the objective function coefficient
         * @return Objective coefficient (negative for biomass reactions)
         */
        double obj_coeff() const {return obj_coeff_;}
        /**
         * @brief Gets the 1-based LP solver ID
         * @return One-based index for LP solver (index + 1)
         */
        int lps_id() const {return (int)index_ + 1;}
        /**
         * @brief Checks if this is a biomass reaction
         * @return true if objective coefficient is negative, false otherwise
         */
        bool is_biomass() const {return is_biomass_;}
        /**
         * @brief Checks if this is an export reaction
         * @return true if all metabolite coefficients are positive, false otherwise
         */
        bool is_export() const {return is_export_;}
        // Is the reaction an EX or DM reaction -
        // -- EX has no reactants
        // -- DM has no products
        /**
         * @brief Checks if this is an exchange (EX) reaction
         * @return true if the reaction has no reactants (substrates), false otherwise
         */
        bool is_ex() const {return in_mets_.size() == 0;}
        /**
         * @brief Checks if this is a demand (DM) reaction
         * @return true if the reaction has no products, false otherwise
         */
        bool is_dm() const {return out_mets_.size() == 0;}
        /**
         * @brief Checks if this is either an exchange or demand reaction
         * @return true if the reaction is EX or DM, false otherwise
         */
        bool is_ex_dm() const {return is_ex() || is_dm();}
        /**
         * @brief Gets the reduced cost from the last LP solution
         * @return Reduced cost value for this reaction's column variable
         */
        double reduced_cost() const {return reduced_cost_;}
        /**
         * @brief Sets the reduced cost from LP solution
         * @param reduced_cost New reduced cost value
         */
        void set_reduced_cost (double reduced_cost) {
            reduced_cost_ = reduced_cost;
        }
        /**
         * @brief Checks if this is a dummy reaction
         * @return true if this is a dummy reaction, false otherwise
         */
        bool is_dummy() const {return is_dummy_;}
        /**
         * @brief Checks if the reaction is active in the network
         * @return true if active, false otherwise
         */
        bool is_active() const {return is_active_;}
        /**
         * @brief Sets the active status of the reaction
         * @param is_active New active status
         */
        void set_active (bool is_active) {
            is_active_ = is_active;
        }
        /**
         * @brief Checks if the reaction is selected for LP optimization
         * @return true if selected, false otherwise
         */
        bool is_selected() const {return is_selected_;}
        /**
         * @brief Sets the selection status for LP optimization
         * @param selected New selection status
         */
        void set_selected (bool selected) {
            DEBUG ('c', "Set select for " << name_ << " to " << selected);
            is_selected_ = selected;
        }
        /**
         * @brief Gets the upper bound on reaction flux
         * @return Flux upper bound value
         */
        double flux_ub() const {return flux_ub_;}
        /**
         * @brief Sets the upper bound on reaction flux
         * @param flux_ub New flux upper bound value
         */
        void set_flux_ub (double flux_ub) {
            flux_ub_ = flux_ub;
        }
        /**
         * @brief Checks if a known flux measurement exists for this reaction
         * @return true if flux is known, false otherwise
         */
        bool has_known_flux() const {return has_known_flux_;}
        /**
         * @brief Gets the known flux measurement
         * @return Known flux value (experimental or reference)
         */
        double known_flux() const {
            //assert (has_known_flux_);
            return known_flux_;
        }
        /**
         * @brief Gets the error/uncertainty in the known flux
         * @return Known flux error value
         */
        double known_flux_error() const {return known_flux_error_; }
        /**
         * @brief Sets the known flux measurement and its error
         * @param known_flux Measured flux value
         * @param known_flux_error Uncertainty in the flux measurement
         */
        void set_known_flux (double known_flux, double known_flux_error) {
            has_known_flux_ = true;
            known_flux_ = known_flux;
            known_flux_error_ = known_flux_error;
        }

        /**
         * @brief Gets all metabolite indices participating in this reaction
         * @return List of metabolite indices
         */
        const std::list<size_t>& mets() const {return mets_;}
        /**
         * @brief Gets input metabolites (substrates/reactants)
         * @return List of metabolite pointers with negative coefficients
         */
        const std::list<const Metabolite*>& in_mets() const {return in_mets_;}
        /**
         * @brief Gets output metabolites (products)
         * @return List of metabolite pointers with positive coefficients
         */
        const std::list<const Metabolite*>& out_mets() const {return out_mets_;}

        /**
         * @brief Sets the stoichiometric coefficient for a metabolite
         *
         * Adds the metabolite to the reaction with the given coefficient.
         * Negative coefficients indicate substrates, positive indicate products.
         * Updates the reaction formula string accordingly.
         *
         * @param met Pointer to the metabolite
         * @param coeff Stoichiometric coefficient (negative for substrates, positive for products)
         */
        void set_coeff (const Metabolite* met, double coeff)
        {
            mets_.push_back (met->index());
            met_coeff_[met->index()] = coeff;
            if (coeff < 0) {
                in_mets_.push_back (met);
                is_export_ = false;
                auto sep = "";
                if (formula_from_.length() > 0)
                    sep = "+ ";
                if (limeDblEqual (coeff, -1.0f))
                    formula_from_ += sep + met->name() + " ";
                else {
                    double rounded = -coeff;
                    if (coeff < -0.01f)
                        limeRound(-coeff * 100.0f) / 100.0f;
                    formula_from_ +=
                        sep + limeFormat("%g", rounded) + " " +
                        met->name() + " ";
                }
            }
            else {
                out_mets_.push_back (met);
                auto sep = " ";
                if (formula_to_.length() > 0)
                    sep = " + ";
                if (limeDblEqual (coeff, 1.0f))
                    formula_to_ += sep + met->name();
                else {
                    double rounded = coeff;
                    if (coeff > 0.01f)
                        limeRound(coeff * 100.0f) / 100.0f;
                    formula_to_ +=
                        sep + limeFormat("%g", rounded) + " " + met->name();
                }
            }
        }
        /**
         * @brief Removes a metabolite from the reaction
         *
         * Removes the metabolite and its coefficient, updates internal lists,
         * and regenerates the formula string.
         *
         * @param met Pointer to the metabolite to remove
         */
        void remove_met (const Metabolite* met)
        {
            mets_.remove (met->index());
            met_coeff_[met->index()] = 0.0f;
            double coeff = met_coeff (met);
            if (coeff < 0.0f) {
                in_mets_.remove (met);
            }
            else if (coeff > 0.0f) {
                out_mets_.remove (met);
            }
            make_formula();
        }
        /**
         * @brief Regenerates the formula string from current metabolites and coefficients
         *
         * Reconstructs the reaction formula by iterating through input and output
         * metabolites with their stoichiometric coefficients.
         */
        void make_formula()
        {
            formula_from_ = "";
            formula_to_ = "";

            for (auto met : in_mets_) {
                double coeff = met_coeff_[met->index()];
                assert (coeff < 0.0f);
                auto sep = "";
                if (formula_from_.length() > 0)
                    sep = "+ ";
                if (limeDblEqual (coeff, -1.0f))
                    formula_from_ += sep + met->name() + " ";
                else {
                    double rounded = -coeff;
                    if (coeff < -0.01f)
                        limeRound(-coeff * 100.0f) / 100.0f;
                    formula_from_ +=
                        sep + limeFormat("%g", rounded) + " " +
                        met->name() + " ";
                }
            }
            for (auto met : out_mets_) {
                double coeff = met_coeff_[met->index()];
                assert (coeff > 0.0f);
                auto sep = " ";
                if (formula_to_.length() > 0)
                    sep = " + ";
                if (limeDblEqual (coeff, 1.0f))
                    formula_to_ += sep + met->name();
                else {
                    double rounded = coeff;
                    if (coeff > 0.01f)
                        limeRound(coeff * 100.0f) / 100.0f;
                    formula_to_ +=
                        sep + limeFormat("%g", rounded) + " " + met->name();
                }
            }
        }
        
        /**
         * @brief Gets the number of metabolites in the reaction
         * @return Count of metabolites with non-zero coefficients
         */
        size_t num_mets() const
        {
            return met_coeff_.size();
        }
        /**
         * @brief Checks if a metabolite has a non-zero coefficient
         * @param k Metabolite index
         * @return true if the metabolite participates in this reaction, false otherwise
         */
        bool met_has_coeff (size_t k) const
        {
            auto iter = met_coeff_.find(k);
            return (iter != met_coeff_.end());
        }
        /**
         * @brief Checks if the reaction uses the given metabolite
         * @param met Pointer to the metabolite
         * @return true if the metabolite has a non-zero coefficient, false otherwise
         */
        bool uses (const Metabolite* met) const
        {
            return met_has_coeff (met->index());
        }
        /**
         * @brief Gets the stoichiometric coefficient for a metabolite by index
         * @param k Metabolite index
         * @return Stoichiometric coefficient (0.0 if metabolite not in reaction)
         */
        double met_coeff (size_t k) const
        {
            auto iter = met_coeff_.find(k);
            if (iter == met_coeff_.end())
                return 0.0f;
            return iter->second;
        }
        /**
         * @brief Gets the stoichiometric coefficient for a metabolite
         * @param met Pointer to the metabolite
         * @return Stoichiometric coefficient (0.0 if metabolite not in reaction)
         */
        double met_coeff (const Metabolite* met) const
        {
            return met_coeff(met->index());
        }
        
        /**
         * @brief Writes the reaction to an output stream with custom cost
         *
         * Outputs the reaction in a specific format including name, cost,
         * flux bounds, metabolites with coefficients, and full name if different.
         *
         * @param out Output stream to write to
         * @param cost Cost value to write (may differ from obj_coeff)
         */
        void write_to (std::ostream& out, double cost) const
        {
            out << "REACTION " << name() <<
                " " << cost <<
                " " << flux_ub();
            for (auto met : in_mets()) {
                out << " " << met->name() <<
                    " " << met_coeff(met);
            }
            for (auto met : out_mets()) {
                out << " " << met->name() <<
                    " " << met_coeff(met);
            }
            if (full_name_.compare (name_) != 0)
                out << " __fullname__ \"" << full_name_ << "\"";
            out << std::endl;
        }

        /**
         * @brief Writes the reaction to an output stream
         *
         * Uses the current objective coefficient as the cost.
         *
         * @param out Output stream to write to
         */
        void write_to (std::ostream& out) const
        {
            write_to (out, obj_coeff());
        }

        /**
         * @brief Checks if this reaction is equivalent to another
         *
         * Compares metabolites and their coefficients to determine if two
         * reactions represent the same biochemical transformation. Comparison
         * is done by metabolite names to allow cross-scenario comparison.
         *
         * @param other Pointer to the reaction to compare with
         * @return true if reactions have identical metabolites and coefficients, false otherwise
         */
        bool same_as (const Reaction* other) const
        {
            if (in_mets_.size() != other->in_mets_.size())
                return false;
            if (out_mets_.size() != other->out_mets_.size())
                return false;
            // Check mets are the same. Note that since these may be
            // different scenarios, we need to check names
            for (auto met1 : in_mets_) {
                bool ok = false;
                for (auto met2 : other->in_mets_) {
                    if (met1->name() == met2->name()) {
                        ok =
                            limeDblEqual (
                                met_coeff (met1), other->met_coeff(met2)
                            );
                        break;
                    }
                }
                if (!ok)
                    return false;
            }
            // Repeat for out-mets
            for (auto met1 : out_mets_) {
                bool ok = false;
                for (auto met2 : other->out_mets_) {
                    if (met1->name() == met2->name()) {
                        ok =
                            limeDblEqual (
                                met_coeff (met1), other->met_coeff(met2)
                            );
                        break;
                    }
                }
                if (!ok)
                    return false;
            }
            return true;
        }

        /**
         * @brief Checks if this reaction is the reverse of another
         *
         * Determines if this reaction represents the reverse direction of
         * another reaction by comparing substrates with products and checking
         * for opposite coefficient signs.
         *
         * @param other Pointer to the reaction to compare with
         * @return true if this reaction is the reverse of the other, false otherwise
         */
        bool reverse_of (const Reaction* other) const
        {
            if (in_mets_.size() != other->out_mets_.size())
                return false;
            if (out_mets_.size() != other->in_mets_.size())
                return false;
            // Check mets are the same. Note that since these may be
            // different scenarios, we need to check names
            DEBUG ('q', name() << " " << other->name() << " in");
            for (auto met1 : in_mets_) {
                bool ok = false;
                for (auto met2 : other->out_mets_) {
                    if (met1->name() == met2->name()) {
                        DEBUG (
                            'q', "  Met " << met1->name() <<
                            " coeff1 " << met_coeff (met1) <<
                            " coeff2 " << -other->met_coeff (met2)
                        );
                        ok =
                            limeDblEqual (
                                met_coeff (met1), 0.0f - other->met_coeff(met2)
                            );
                        break;
                    }
                }
                DEBUG ('q', "    OK " << ok);
                if (!ok)
                    return false;
            }
            // Repeat for out-mets
            DEBUG ('q', name() << " " << other->name() << " out");
            for (auto met1 : out_mets_) {
                bool ok = false;
                for (auto met2 : other->in_mets_) {
                    if (met1->name() == met2->name()) {
                        DEBUG (
                            'q', "  Met " << met1->name() <<
                            " coeff1 " << met_coeff (met1) <<
                            " coeff2 " << -other->met_coeff (met2)
                        );
                        ok =
                            limeDblEqual (
                                met_coeff (met1), 0.0f - other->met_coeff(met2)
                            );
                        break;
                    }
                }
                DEBUG ('q', "    OK " << ok);
                if (!ok)
                    return false;
            }
            return true;
        }

    
    private:
        std::string name_; /**< Short name identifier for the reaction */
        std::string full_name_; /**< Full descriptive name of the reaction */
        std::string formula_from_; /**< String representation of substrates (left side of reaction) */
        std::string formula_to_; /**< String representation of products (right side of reaction) */
        size_t index_; /**< Zero-based index of this reaction in the network */
        double obj_coeff_; /**< Objective function coefficient (negative for biomass reactions) */
        bool is_biomass_; /**< True if this is a biomass reaction (obj_coeff < 0) */
        bool is_export_; /**< True if this is an export reaction (all coefficients positive) */
        std::list<size_t> mets_; /**< Indices of metabolites participating in this reaction */
        std::list<const Metabolite*> in_mets_; /**< Substrates/reactants (negative coefficients) */
        std::list<const Metabolite*> out_mets_; /**< Products (positive coefficients) */
        std::map<size_t,double> met_coeff_; /**< Map of metabolite index to stoichiometric coefficient */
        double flux_ub_; /**< Upper bound on reaction flux */
        bool has_known_flux_; /**< True if a known flux measurement exists */
        double known_flux_; /**< Known/measured flux value (negative if unknown) */
        double known_flux_error_; /**< Uncertainty/error in the known flux measurement */ 

        bool is_dummy_; /**< True if this is a dummy reaction for a particular metabolite */
        bool is_active_; /**< True if this reaction is active in the network */
        bool is_selected_; /**< True if this reaction is selected for LP optimization */
        double reduced_cost_; /**< Reduced cost of the column representing this reaction in the last LP */
    };
    
    /**
     * @typedef ReactionPtr
     * @brief Shared pointer type for Reaction objects
     */
    using ReactionPtr = std::shared_ptr<Reaction>;
}
