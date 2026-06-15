/**
 * @file metabolite.h
 * @brief Defines the Metabolite class for representing biochemical compounds in metabolic networks.
 *
 * This file contains the Metabolite class used in metabolic flux balance analysis and gap-filling.
 * Metabolites represent chemical species that participate in biochemical reactions within
 * metabolic networks.
 */
#pragma once

#include <iostream>
#include <memory>

#include "lime/numutil.h"
#include "lime/displayable.h"

/**
 * @brief Namespace for metabolic optimization and simulation harness components.
 */
namespace mosh
{
    /**
     * @class Metabolite
     * @brief Represents a metabolite (biochemical compound) in a metabolic network.
     *
     * The Metabolite class encapsulates information about chemical species in metabolic
     * pathways, including identification, categorization, and special properties relevant
     * to flux balance analysis and gap-filling algorithms. Metabolites can be designated
     * as sources, carbon sources, cycle metabolites, or dummy metabolites used for
     * computational purposes.
     */
    class Metabolite : public lime::Displayable
    {
    public:
        /**
         * @brief Constructs a new Metabolite object.
         * @param name Short identifier for the metabolite (e.g., "glc_D" for D-glucose).
         * @param full_name Full descriptive name of the metabolite (e.g., "D-Glucose").
         */
        Metabolite(std::string name, std::string full_name) :
            Displayable(),
            name_(name),
            full_name_(full_name),
            is_source_(false),
            is_cycle_met_(false),
            is_c_source_(false),
            is_dummy_(false)
        {
        }

        /**
         * @brief Gets the short identifier name of the metabolite.
         * @return The metabolite's short name.
         */
        std::string name() const override {return name_;}
        /**
         * @brief Sets the short identifier name of the metabolite.
         * @param name The new short name for the metabolite.
         */
        void set_name(std::string name) {name_ = name;}
        /**
         * @brief Gets the full descriptive name of the metabolite.
         * @return The metabolite's full name.
         */
        std::string full_name() const {return full_name_;}
        /**
         * @brief Gets the index of this metabolite in the stoichiometric matrix.
         * @return The zero-based index position.
         */
        size_t index() const {return index_;}
        /**
         * @brief Sets the index of this metabolite in the stoichiometric matrix.
         * @param index The zero-based index position to assign.
         */
        void set_index (size_t index) {index_ = index;}
        /**
         * @brief Checks if this metabolite is designated as a source compound.
         * @return True if the metabolite is a source, false otherwise.
         */
        bool is_source() const {return is_source_;}
        /**
         * @brief Sets whether this metabolite is a source compound.
         * @param is_source True to designate as a source, false otherwise.
         */
        void set_is_source (bool is_source) {
            is_source_ = is_source;
        }
        /**
         * @brief Checks if this metabolite participates in metabolic cycles.
         * @return True if the metabolite is part of a cycle, false otherwise.
         */
        bool is_cycle_met() const {return is_cycle_met_;}
        /**
         * @brief Sets whether this metabolite participates in metabolic cycles.
         * @param is_cycle_met True if part of a cycle, false otherwise.
         */
        void set_cycle_met (bool is_cycle_met) {
            is_cycle_met_ = is_cycle_met;
        }
        /**
         * @brief Checks if this metabolite is a carbon source.
         * @return True if the metabolite is a carbon source, false otherwise.
         */
        bool is_c_source() const {return is_c_source_;}
        /**
         * @brief Sets whether this metabolite is a carbon source.
         * @param is_c_source True to designate as a carbon source, false otherwise.
         */
        void set_c_source (bool is_c_source) {
            is_c_source_ = is_c_source;
        }
        /**
         * @brief Checks if this is a dummy metabolite used for computational purposes.
         * @return True if the metabolite is a dummy, false otherwise.
         */
        bool is_dummy() const {return is_dummy_;}
        /**
         * @brief Sets whether this is a dummy metabolite.
         * @param is_dummy True to designate as a dummy metabolite, false otherwise.
         */
        void set_is_dummy (bool is_dummy) {
            is_dummy_ = is_dummy;
        }
        
        /**
         * @brief Writes the metabolite information to an output stream.
         * @param out The output stream to write to.
         */
        void write_to (std::ostream& out) const
        {
            out << "MET " << name() <<
                " \"" << full_name() << "\"" << std::endl;
        }
            
    private:
        std::string name_; /**< Short identifier name of the metabolite */"
        std::string full_name_; /**< Full descriptive name of the metabolite */"
        size_t index_; /**< Index position in the stoichiometric matrix */"
        bool is_source_; /**< Flag indicating if this is a source metabolite */"
        bool is_cycle_met_; /**< Flag indicating if this metabolite participates in cycles */"
        bool is_c_source_; /**< Flag indicating if this is a carbon source */"
        bool is_dummy_; /**< Flag indicating if this is a dummy metabolite for computational purposes */"
    };
    
    /**
     * @typedef MetabolitePtr
     * @brief Shared pointer type for Metabolite objects.
     */
    using MetabolitePtr = std::shared_ptr<Metabolite>;
}
