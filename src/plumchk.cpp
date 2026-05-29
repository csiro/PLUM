#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <set>

#include <stdlib.h>

#include "lime/opts.h"
#include "lime/dig.h"
#include "lime/debug.h"
#include "lime/error.h"
#include "lime/strutil.h"
#include "lime/lockfile.h"

#include "mosh/builddate.h" // For _build_date
#include "mosh/scenario.h"
#include "mosh/params.h"

using namespace std;
using namespace lime;
using namespace mosh;

// Define the OS variable
const char* 
get_os_string() {
#ifdef __CYGWIN__
    return "Cygwin";
#elif _WIN32
    return "Windows";
#elif __APPLE__
    return "Mac OS X";
#elif __linux__
    return "Linux";
#elif __unix__
    return "Unix";
#else
    return "Unknown OS";
#endif
}

#define DISPLAY(X) {if (verbose) {cout << X << endl;} if (out != nullptr) { *out << X << endl;}}

int
main (int argc, const char* argv[]) 
{
    int seed = 0;
    string scenario1_fn = "";
    string scenario2_fn = "";
    string merge_fn = "";
    string out_fn = "";
    string summary_fn = "";
    string debug_str = "";
    bool keep_rev = false;
    bool verbose = true;

    Opts opts (R"EOF(
    Identify duplicate or reversed reactions.
    Can create merged file that has all reactions in scenario 1, plus
    those reactions in scenario 2 that
    - Do not appear in scenario 1
    - Are not the reverse of a reaction in scenario 1
    (Logic: If you wanted the reverse in scenario 1, then you would have put it in)
)EOF"
    );

    opts.add_opt ("-q", &verbose, "Quiet (don't echo to stdout)");
    opts.add_opt_filename ("-o", &out_fn, "Output filename");
    opts.add_opt_filename ("-m", &merge_fn, "Merged output filename");
    opts.add_opt_filename ("-O", &summary_fn, "Summary filename (appended)");
    opts.add_opt ("-rev", &keep_rev, "Keep reverse reactions in scenario 2");
    opts.add_opt ("-d", &debug_str, "Debug string");

    // Add a filename args 
    opts.add_arg ("scenario1.dat", &scenario1_fn, "Scenario 1");
    opts.add_arg ("scenario2.dat", &scenario2_fn, "Scenario 2");
    
    string debug = "DEBUG";
#ifdef NDEBUG
    debug = "Release";
#endif
    const char* os = get_os_string();
    stringstream prog_stream;
    prog_stream << argv[0] << " Version of " << _build_date << " " <<
        buildId() << " " << os << " " << debug;
    string progname = prog_stream.str();
    cerr << progname << endl;

    if (!opts.process (argc, argv))
        exit (1);
    
#ifdef NDEBUG
    if (debug_str.compare("") != 0)
        cerr << "Debug string ignored in Release version" << endl;
#else
    Debug::setFilename ("debug.out");
    Debug::setKey (debug_str);
    Debug::debugFile() << progname << endl;
    Debug::debugFile() << "Timestamp " << lime::todayString() << endl;
    Debug::debugFile() << "DebugStr: " << debug_str << endl;
#endif

    shared_ptr<ofstream> out = nullptr;
    if (out_fn.length() > 0) {
        cout << "Writing " << out_fn << endl;
        out = make_shared<ofstream> (out_fn);
        *out << "# Produced by " << progname << " on " << todayString() << endl;
        *out << "# scenario1 "  << scenario1_fn <<
            " scenario2 "  << scenario2_fn << endl;
    }

    Params params;
    
    Scenario scenario1;
    Scenario scenario2;

    scenario1.read_data (scenario1_fn);
    scenario2.read_data (scenario2_fn);

    for (auto& react2 : scenario2.reactions()) {
        react2->set_selected(true);
    }
    
    using ReactPair = pair<const Reaction*, const Reaction*>;
    list<ReactPair> scen1_repeats;	// Repeats, scenario 1
    list<ReactPair> scen2_repeats;	// Repeats, scenario 2
    list<ReactPair> same_diff; 	// Same name different formula
    list<ReactPair> diff_same; 	// Different name same formula
    list<ReactPair> rev;	// reverse
    list<ReactPair> unex_rev;	// unexpected reverse
    size_t same_count = 0;
    size_t dropped = 0;

    // Check for and remove internal repeats in scenario 1
    for (size_t i = 0; i < scenario1.num_reactions(); i++) {
        auto react_i = scenario1.reaction(i);
        if (!react_i->is_selected()) // Already dropped
            continue; 
        DEBUG ('c', "React1_i " << *react_i);
        for (size_t j = i+1; j < scenario1.num_reactions(); j++) {
            auto react_j = scenario1.reaction(j);
            if (!react_j->is_selected()) // Already dropped
                continue; 
            if (react_i->same_as (react_j)) {
                // Save cheaper version
                if (react_i->obj_coeff() <= react_j->obj_coeff()) {
                    DEBUG (
                        'c', "    React1_j " << *react_j <<
                        " " << react_j->obj_coeff() << 
                        " is more expensive repeat of React1_i " << *react_i <<
                        " " << react_i->obj_coeff() 
                    );
                    react_j->set_selected (false);
                }
                else {
                    DEBUG (
                        'c', "    React1_j " << *react_j <<
                        " " << react_j->obj_coeff() << 
                        " is cheaper repeat of React1_i " << *react_i <<
                        " " << react_i->obj_coeff() 
                    );
                    react_i->set_selected (false);
                }
                scen1_repeats.push_back (make_pair(react_i,react_j));
                dropped++;
            }
        }
    }
    
    // Check for and remove internal repeats in scenario 2
    for (size_t i = 0; i < scenario2.num_reactions(); i++) {
        auto react_i = scenario2.reaction(i);
        if (!react_i->is_selected()) // Already dropped
            continue; 
        DEBUG ('c', "React2_i " << *react_i);
        for (size_t j = i+1; j < scenario2.num_reactions(); j++) {
            auto react_j = scenario2.reaction(j);
            if (!react_j->is_selected()) // Already dropped
                continue; 
            if (react_i->same_as (react_j)) {
                // Save cheaper version
                if (react_i->obj_coeff() <= react_j->obj_coeff()) {
                    DEBUG (
                        'c', "    React1_j " << *react_j <<
                        " " << react_j->obj_coeff() << 
                        " is more expensive repeat of React1_i " << *react_i <<
                        " " << react_i->obj_coeff() 
                    );
                    react_j->set_selected (false);
                }
                else {
                    DEBUG (
                        'c', "    React1_j " << *react_j <<
                        " " << react_j->obj_coeff() << 
                        " is cheaper repeat of React1_i " << *react_i <<
                        " " << react_i->obj_coeff() 
                    );
                    react_i->set_selected (false);
                }
                scen2_repeats.push_back (make_pair(react_i, react_j));
                dropped++;
            }
        }
    }
    
    for (auto& react1 : scenario1.reactions()) {
        DEBUG ('c', "React1 " << *react1);
        if (!react1->is_selected()) // Already dropped
            continue; 
        for (auto& react2 : scenario2.reactions()) {
            if (!react2->is_selected()) // Already dropped
                continue; 
            bool same_name = (react1->name() == react2->name());
            bool same_formula = react1->same_as (react2.get());
            bool is_rev = react1->reverse_of (react2.get());

            if (same_formula)  {
                react2->set_selected (false);
                same_count++;
                dropped++;
            }
            if (same_name && !same_formula) {
                same_diff.push_back (make_pair(react1.get(), react2.get()));
                react2->set_selected (false);
                dropped++;
            }
            if (!same_name && same_formula) 
                diff_same.push_back (make_pair(react1.get(), react2.get()));
            
            DEBUG (
                'c', "  React2 " << *react2 <<
                " same_name " << same_name <<
                " same_formula " << same_formula <<
                " is_rev " << is_rev <<
                " is_sel " << react2->is_selected()
            );
            if (is_rev && react2->is_selected()) {
                if (!keep_rev) {
                    rev.push_back (make_pair(react1.get(), react2.get()));
                    react2->set_selected (false);
                    dropped++;
                }

                if (
                    (react1->name() != (react2->name() + "_rev")) &&
                    (react2->name() != (react1->name() + "_rev"))
                ) {
                    // Make sure it's not an ex/dm pair
                    bool exdm1 = 
                        react1->name().compare(0, 3, string("EX_")) == 0 ||
                        react1->name().compare(0, 3, string("DM_")) == 0;
                    bool exdm2 = 
                        react2->name().compare(0, 3, string("EX_")) == 0 ||
                        react2->name().compare(0, 3, string("DM_")) == 0;
                    bool rest =
                        react1->name().substr(2).compare(
                            react2->name().substr(2)
                        ) == 0;
                    if (!exdm1 || !exdm2 || !rest) {
                        unex_rev.push_back (
                            make_pair(react1.get(), react2.get())
                        );
                    }
                }
            }
        }
    }
    DISPLAY ("Scenario 1 repeats");
    for (auto& p : scen1_repeats) {
        auto react1 = p.first;
        auto react2 = p.second;
        DISPLAY ("  " << react1->name() << " " << react2->name());
    }
    DISPLAY ("");
    DISPLAY ("Scenario 2 repeats");
    for (auto& p : scen2_repeats) {
        auto react1 = p.first;
        auto react2 = p.second;
        DISPLAY ("  " << react1->name() << " " << react2->name());
    }
    DISPLAY ("");
    DISPLAY ("Same name different formula");
    for (auto& p : same_diff) {
        auto react1 = p.first;
        auto react2 = p.second;
        DISPLAY ("  " << react1->name() << " " << react2->name());
    }
    DISPLAY ("");
    DISPLAY ("Different name same formula");
    for (auto& p : diff_same) {
        auto react1 = p.first;
        auto react2 = p.second;
        DISPLAY ("  " << react1->name() << " " << react2->name());
    }
    DISPLAY ("");
    DISPLAY ("Unexpected Reverse");
    for (auto& p : unex_rev) {
        auto react1 = p.first;
        auto react2 = p.second;
        DISPLAY ("  " << react1->name() << " " << react2->name());
    }
    DISPLAY ("");
    DISPLAY ("Reverse");
    for (auto& p : rev) {
        auto react1 = p.first;
        auto react2 = p.second;
        DISPLAY ("  " << react1->name() << " " << react2->name());
    }
    
    stringstream summary;
    summary <<
        " scenario1 " << scenario1_fn <<
        " mets1 " << scenario1.num_metabolites() <<
        " reactions1 " << scenario1.num_reactions() <<
        " scenario2 " << scenario2_fn <<
        " mets2 " << scenario2.num_metabolites() <<
        " reactions2 " << scenario2.num_reactions() <<
        " scen1_repeats " << scen1_repeats.size() <<
        " scen2_repeats " << scen2_repeats.size() <<
        " same_count " << same_count <<
        " rev_count " << rev.size() <<
        " dropped " << dropped << 
        " same_diff_count " << same_diff.size() <<
        " diff_same_count " << diff_same.size() <<
        " build " << buildId();

    cout << summary.str() << endl;

    if (merge_fn.length() > 0) {
        cout << "Writing " << merge_fn << endl;
        ofstream out (merge_fn);
        out << "# Produced by " << progname << " on " << todayString() << endl;
        out << "# By merging "  << scenario1_fn <<
            " and "  << scenario2_fn << endl;
        out << "# " << summary.str() << endl;
        // Write scen 1 mets
        for (auto& met : scenario1.metabolites()) {
            met->write_to(out);
        }
        vector<bool> seen (scenario2.num_metabolites(), false);
        // Add unseen scen 2 mets
        for (auto& react2 : scenario2.reactions()) {
            if (react2->is_selected()) {
                for (auto k : react2->mets()) {
                    if (seen[k])
                        continue;
                    seen[k] = true; // Don't have to process this met again
                    auto met = scenario2.metabolite(k);
                    if (scenario1.find_metabolite (met->name()) == nullptr) {
                        // Add this met
                        met->write_to (out);
                    }
                }
            }
        }
        // Write scen 1 reactions
        for (auto& react1 : scenario1.reactions()) {
            if (react1->is_selected()) {
                react1->write_to (out);
            }
        }
        // Write selected scen 2 reactions
        for (auto& react2 : scenario2.reactions()) {
            if (react2->is_selected()) {
                react2->write_to (out);
            }
        }
    }
    
    if (out_fn.length() > 0) {
        cout << "Wrote " << out_fn << endl;
        out->close();
    }
    if (summary_fn.length() > 0) {
        cout << "Writing " << summary_fn << endl;
        LockFile lock (summary_fn);
        ofstream out (summary_fn, fstream::out | fstream::app);
        out << summary.str() << endl;
    }

    return 0;
}
