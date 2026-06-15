/**
 * @file mergeplum.cpp
 * @brief Merges two PLUM metabolic network scenarios into a single scenario
 *
 * This utility merges two scenario files (.pld format) used in metabolic gap-filling
 * analysis. It handles the reconciliation of metabolites and reactions between scenarios,
 * validating stoichiometry consistency and allowing selective preservation of reaction
 * costs from either source scenario. The merged output maintains metabolic network
 * integrity while combining unique elements from both inputs.
 */
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>

#include <stdlib.h>

#include "lime/opts.h"
#include "lime/dig.h"
#include "lime/debug.h"
#include "lime/linereader.h"
#include "lime/limetok.h"
#include "lime/strutil.h"
#include "lime/sortpairt.h"

#include "mosh/builddate.h" // For _build_date
#include "mosh/gapsolver.h"
#include "mosh/params.h"

using namespace std;
using namespace lime;
using namespace mosh;

// Define the OS variable
/**
 * @brief Detects and returns the operating system platform at compile time
 *
 * Uses preprocessor directives to identify the compilation platform and return
 * a human-readable string representation.
 *
 * @return Constant string identifying the OS: "Cygwin", "Windows", "Mac OS X",
 *         "Linux", "Unix", or "Unknown OS"
 */
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

/**
 * @brief Main entry point for the mergeplum scenario merging utility
 *
 * Merges two PLUM scenario files by combining their metabolites and reactions.
 * The merge process:
 * - Validates that reactions with identical names have matching stoichiometry
 * - Identifies metabolites that appear in both scenarios using name matching
 * - Uses the -w flag to determine which scenario's scores (reaction costs) to preserve
 * - Writes unused metabolites and unique reactions from each scenario
 * - Generates warnings for metabolite full name mismatches
 * - Generates errors for stoichiometry inconsistencies
 * - Outputs a merged scenario file and a diagnostic report
 *
 * @param argc Number of command-line arguments
 * @param argv Array of command-line argument strings
 * @return 0 on successful merge, 1 on command-line parsing failure
 */
int
main (int argc, const char* argv[]) 
{
    string scen1_fn = "";
    string scen2_fn = "";
    string out_fn = "out.pld";
    string debug_str = "";

    // Keep scores from which scenario
    int which_scen = 1;

    Params params;

    Opts opts (R"EOF(
    Merge two plum scenarios.
    - Reactions with the same name are checked for equality
    - If both scenarios contain a metabolite or reaction,
      '-w' switch controls which version is used
    - In particular '-w' switch controls reaction costs.
)EOF"
    );

    opts.add_opt ("-w", &which_scen, "Use scores from which scenario");
    opts.add_opt ("-d", &debug_str, "Debug string");

    opts.add_arg ("scen1.pld", &scen1_fn, "Scenario 1");
    opts.add_arg ("scen2.pld", &scen2_fn, "Scenario 2");
    opts.add_optional_arg ("out.pld", &out_fn, "Merged scenario");
    
    const char* debug = "DEBUG";
#ifdef NDEBUG
    debug = "Release";
#endif
    const char* os = get_os_string();
    stringstream prog_stream;
    prog_stream << argv[0] << " Version of " << _build_date << " " <<
        os << " " << debug;
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

    if (which_scen != 1 && which_scen != 2)
        opts.usage ("-w must be 1 or 2");

    Scenario scen1;
    Scenario scen2;
    
    cout << "Reading data from " << scen1_fn << endl;
    scen1.read_data (scen1_fn);

    cout << "Reading data from " << scen2_fn << endl;
    scen2.read_data (scen2_fn);

    cout << "Compare" << endl;

    list<string> warn;
    int warn_count = 0;
    list<string> error;
    int error_count = 0;

    // Translation from scen1 metabolite to scen2
    vector<const Metabolite*> xlate_met_12 (scen1.num_metabolites(), nullptr);
    // Translation from scen2 metabolite to scen1
    vector<const Metabolite*> xlate_met_21 (scen2.num_metabolites(), nullptr);

    for (auto& met1 : scen1.metabolites()) {
        auto met2 = scen2.find_metabolite (met1->name());
        xlate_met_12[met1->index()] = met2;
        if (met2 != nullptr)
            xlate_met_21[met2->index()] = met1.get();
        if (met2 != nullptr) {
            size_t len1 = met1->full_name().length();
            size_t len2 = met2->full_name().length();
            sizet minlen = std::min(len1,len2);
            if (
                !equal_ic (
                    met1->full_name().substr(0,minlen),
                    met2->full_name().substr(0,minlen)
                )
            ) {
                warn_count++;
                warn.push_back (
                    "Metabolite " + met1->name() + " differs in full name "
                );
                warn.push_back (
                    limeFormat (
                        "%20s: %s", scen1_fn.c_str(), met1->full_name().c_str()
                    )
                );
                warn.push_back (
                    limeFormat (
                        "%20s: %s", scen2_fn.c_str(), met2->full_name().c_str()
                    )
                warn.push_back (
                    limeFormat (
                        "%20s: %s", scen1_fn.c_str(), met1->full_name().substr(0,minlen).c_str()
                    )
                );
                warn.push_back (
                    limeFormat (
                        "%20s: %s", scen2_fn.c_str(), met2->full_name().substr(0,minlen).c_str()
                    )
                );
            }
        }
    }

    vector<bool> met_used1 (scen1.num_metabolites(), false);
    vector<bool> met_used2 (scen2.num_metabolites(), false);

    for (auto& react1 : scen1.reactions()) {
        auto react2 = scen2.find_reaction (react1->name());

        if (react2 == nullptr)
            continue; // No overlap
        
        bool ok = true;
        string error_met = "";
        // Check all react1 mets appear in react2, and have the same coeff
        for (auto idx1 : react1->mets()) {
            auto met1 = scen1.metabolite(idx1);
            met_used1[idx1] = true;

            auto met2 = xlate_met_12[idx1];
            if (met2 == nullptr) {
                ok = false;
                error_met = met1->name();
                break;
            }
            else {
                double coeff1 = react1->met_coeff (met1);
                double coeff2 = react2->met_coeff (met2);
                if (!limeDblEqual (coeff1, coeff2)) {
                    ok = false;
                    error_met = met1->name();
                    break;
                }
            }
        }
        if (ok) {
            // Check all react2 mets appear in react1
            // (already checked the coeff if they both appear)
            for (auto idx2 : react2->mets()) {
                auto met2 = scen2.metabolite(idx2);
                
                auto met1 = xlate_met_21[idx2];
                if (met1 == nullptr) {
                    ok = false;
                    error_met = met2->name();
                    break;
                }
                else if (!react1->uses (met1)) {
                    ok = false;
                    error_met = met1->name();
                    break;
                }
            }
        }
        if (!ok) {
            error_count++;
            error.push_back (
                "Reaction " + react1->name() +
                " differs in stoichiometry for metabolite " + error_met
            );
            error.push_back (
                limeFormat (
                    "%20s: %s", scen1_fn.c_str(), react1->formula().c_str()
                )
            );
            error.push_back (
                limeFormat (
                    "%20s: %s", scen2_fn.c_str(), react2->formula().c_str()
                )
            );
        }
    }

    for (auto& react2 : scen2.reactions()) {
        for (auto idx2 : react2->mets()) {
            met_used2[idx2] = true;
        }
    }
    
    // Write out merged scenario

    ofstream out (out_fn);
    out << "# Produced from " << scen1_fn << " and " <<  scen2_fn <<
        " on " << todayString() << " by " << progname << endl;
    
    for (auto& met1 : scen1.metabolites()) {
        if (!met_used1[met1->index()])
            continue;
        
        if (
            xlate_met_12[met1->index()] == nullptr ||
            which_scen == 1
        ) {
            met1->write_to (out);
        }
        else {
            xlate_met_12[met1->index()]->write_to (out);
        }
    }
    // Write out mets in scen2 that have no equiv in scen1
    for (auto& met2 : scen2.metabolites()) {
        if (!met_used2[met2->index()])
            continue;
        auto met1 = xlate_met_21[met2->index()];
        if (met1 == nullptr)
            met2->write_to (out);
        else {
            // Double check we were not relying on it
            if (!met_used1[met1->index()]) {
                // Huh. Scen 2 uses it but not scen 1. Write it out
                met2->write_to (out);
            }
        }
    }

    for (auto& react1 : scen1.reactions()) {
        auto react2 = scen2.find_reaction (react1->name());
        
        if (react2 == nullptr || which_scen == 1) {
            react1->write_to (out);
        }
        else {
            react2->write_to (out);
        }
    }
    // WRite out reacts from scen2 that don't appear in scen1
    for (auto& react2 : scen2.reactions()) {
        auto react1 = scen1.find_reaction (react2->name());
        
        if (react1 == nullptr) {
            react2->write_to (out);
        }
    }

    string err_fn = "mergeplum.out";
    ofstream errout (err_fn);
    errout << "Produced on " << todayString() << " by " << progname << endl;
        
    cout << "Merge complete" << endl;
    if (error.size() > 0) {
        cout << error_count << " errors" << endl;
        errout << error_count << " errors" << endl;
        for (auto text : error) {
            cout <<  text << endl;
            errout <<  text << endl;
        }
    }
    if (warn.size() > 0) {
        cout << warn_count << " warnings" << endl;
        errout << warn_count << " warnings" << endl;
        for (auto text : warn) {
            cout <<  text << endl;
            errout <<  text << endl;
        }
    }
    cout << endl;
    cout << "Wrote merged scenario to " << out_fn << endl;
    cout << "         Wrote errors to " << err_fn << endl;
    
    return 0;
}
