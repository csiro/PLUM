#pragma once

#include <iostream>
#include <memory>

#include "lime/numutil.h"
#include "lime/displayable.h"

namespace mosh
{
    class Metabolite : public lime::Displayable
    {
    public:
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

        std::string name() const override {return name_;}
        void set_name(std::string name) {name_ = name;}
        std::string full_name() const {return full_name_;}
        size_t index() const {return index_;}
        void set_index (size_t index) {index_ = index;}
        bool is_source() const {return is_source_;}
        void set_is_source (bool is_source) {
            is_source_ = is_source;
        }
        bool is_cycle_met() const {return is_cycle_met_;}
        void set_cycle_met (bool is_cycle_met) {
            is_cycle_met_ = is_cycle_met;
        }
        bool is_c_source() const {return is_c_source_;}
        void set_c_source (bool is_c_source) {
            is_c_source_ = is_c_source;
        }
        bool is_dummy() const {return is_dummy_;}
        void set_is_dummy (bool is_dummy) {
            is_dummy_ = is_dummy;
        }
        
        void write_to (std::ostream& out) const
        {
            out << "MET " << name() <<
                " \"" << full_name() << "\"" << std::endl;
        }
            
    private:
        std::string name_;
        std::string full_name_;
        size_t index_;
        bool is_source_;
        bool is_cycle_met_;
        bool is_c_source_;
        bool is_dummy_;
    };
    
    using MetabolitePtr = std::shared_ptr<Metabolite>;
}
