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

namespace mosh
{
    class Reaction : public lime::Displayable
    {
    public:
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

        std::string name() const override {return name_;}
        void set_name(std::string name) {name_ = name;}
        std::string full_name() const {return full_name_;}
        void set_full_name(std::string full_name) {full_name_ = full_name;}
        std::string formula() const
        {
            return formula_from_ + "->" + formula_to_;
        }
        std::string nice_formula() const
        {
            return nicify (formula_from_) + "->" + nicify (formula_to_);
        }
        std::string nicify (std::string name) const
        {
            std::string nice_name = name;
            lime::replaceAll (nice_name, "__40__", "(");
            lime::replaceAll (nice_name, "__41__", ")");
            lime::replaceAll (nice_name, "__91__", "[");
            lime::replaceAll (nice_name, "__93__", "]");
            return nice_name;
        }
        size_t index() const {return index_;}
        void set_index (size_t index) {index_ = index;}
        void set_obj_coeff (double obj_coeff) {obj_coeff_ = obj_coeff;}
        double obj_coeff() const {return obj_coeff_;}
        int lps_id() const {return (int)index_ + 1;}
        bool is_biomass() const {return is_biomass_;}
        bool is_export() const {return is_export_;}
        // Is the reaction an EX or DM reaction -
        // -- EX has no reactants
        // -- DM has no products
        bool is_ex() const {return in_mets_.size() == 0;}
        bool is_dm() const {return out_mets_.size() == 0;}
        bool is_ex_dm() const {return is_ex() || is_dm();}
        double reduced_cost() const {return reduced_cost_;}
        void set_reduced_cost (double reduced_cost) {
            reduced_cost_ = reduced_cost;
        }
        bool is_dummy() const {return is_dummy_;}
        bool is_active() const {return is_active_;}
        void set_active (bool is_active) {
            is_active_ = is_active;
        }
        bool is_selected() const {return is_selected_;}
        void set_selected (bool selected) {
            DEBUG ('c', "Set select for " << name_ << " to " << selected);
            is_selected_ = selected;
        }
        double flux_ub() const {return flux_ub_;}
        void set_flux_ub (double flux_ub) {
            flux_ub_ = flux_ub;
        }
        bool has_known_flux() const {return has_known_flux_;}
        double known_flux() const {
            //assert (has_known_flux_);
            return known_flux_;
        }
        double known_flux_error() const {return known_flux_error_; }
        void set_known_flux (double known_flux, double known_flux_error) {
            has_known_flux_ = true;
            known_flux_ = known_flux;
            known_flux_error_ = known_flux_error;
        }

        const std::list<size_t>& mets() const {return mets_;}
        const std::list<const Metabolite*>& in_mets() const {return in_mets_;}
        const std::list<const Metabolite*>& out_mets() const {return out_mets_;}

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
        
        size_t num_mets() const
        {
            return met_coeff_.size();
        }
        // Does the met have a non-zero coeff in this reaction?
        bool met_has_coeff (size_t k) const
        {
            auto iter = met_coeff_.find(k);
            return (iter != met_coeff_.end());
        }
        // Does the met have a non-zero coeff in this reaction?
        bool uses (const Metabolite* met) const
        {
            return met_has_coeff (met->index());
        }
        double met_coeff (size_t k) const
        {
            auto iter = met_coeff_.find(k);
            if (iter == met_coeff_.end())
                return 0.0f;
            return iter->second;
        }
        double met_coeff (const Metabolite* met) const
        {
            return met_coeff(met->index());
        }
        
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

        void write_to (std::ostream& out) const
        {
            write_to (out, obj_coeff());
        }

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
        std::string name_;
        std::string full_name_;
        std::string formula_from_;
        std::string formula_to_;
        size_t index_;
        double obj_coeff_;
        // Is this a biomass reaction? true iff obj_coeff < 0
        bool is_biomass_;
        // Is this an export reaction? (all met coeffs positive)
        bool is_export_;
         // The metabolites used in this reaction
        std::list<size_t> mets_;
        std::list<const Metabolite*> in_mets_;
        std::list<const Metabolite*> out_mets_;
        std::map<size_t,double> met_coeff_;
        double flux_ub_;
        bool has_known_flux_;
        double known_flux_; // < 0 means flux not known
        double known_flux_error_; 

        /// Am I a dummy reaction (for a particualr metabolite)
        bool is_dummy_;
        // Am I active 
        bool is_active_;
        // Am I selected for the LP?
        bool is_selected_;
        // Reduced of the col representing this reaction in the last LP
        double reduced_cost_;
    };
    
    using ReactionPtr = std::shared_ptr<Reaction>;
}
