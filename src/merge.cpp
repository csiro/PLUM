/**
 * @file merge.cpp
 * @brief Tool for merging two metabolic gap-filling scenarios
 *
 * This utility merges two scenario files (.pld) by comparing and combining
 * reactions and metabolites. Reactions with the same name are checked for
 * equality, and the user can select which scenario's costs/scores are used
 * when a reaction appears in both scenarios.
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

/**
 * @brief Returns a string identifying the operating system at compile time
 *
 * Uses preprocessor macros to determine the target operating system
 * and returns a human-readable string representation.
 *
 * @return Constant string describing the operating system ("Cygwin", "Windows",
 *         "Mac OS X", "Linux", "Unix", or "Unknown OS")
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
 * @brief Main entry point for the scenario merge tool
 *
 * Reads two metabolic gap-filling scenario files, compares their reactions
 * and metabolites for consistency, and merges them into a single output file.
 * The tool creates bidirectional translation maps between metabolites in the
 * two scenarios and validates reaction equivalence.
 *
 * @param argc Number of command line arguments
 * @param argv Array of command line argument strings
 *             Expected: scen1.pld scen2.pld [out.pld]
 * @return 0 on success, 1 on failure
 */
int
main (int argc, const char* argv[])
{
    /** @brief Filename for the first scenario input file */
    string scen1_fn = "";
    /** @brief Filename for the second scenario input file */
    string scen2_fn = "";
    /** @brief Filename for the merged scenario output file (default: out.pld) */
    string out_fn = "out.pld";
    /** @brief Debug configuration string for diagnostic output */
    string debug_str = "";

    /** @brief Specifies which scenario's scores/costs to use (1 or 2) when reactions appear in both */
    int which_scen = 1;

    /** @brief Parameter configuration object for gap-filling settings */
    Params params;

    /** @brief Command-line options parser for merge tool arguments */
    Opts opts (R"EOF(
    Merge two scenarios.
    - Reactions with the same name are checked for equality
    - Select which scenario's costs are used when a reaction appears in both
)EOF"
    );

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

    assert (data_fn.length() > 0); // Not optional in Opts

    /** @brief First scenario object containing metabolites, reactions, and constraints */
    Scenario scen1;
    /** @brief Second scenario object containing metabolites, reactions, and constraints */
    Scenario scen2;
    
    cout << "Reading data from " << scen1_fn << endl;
    scen1.read_data (scen1_fn);

    cout << "Reading data from " << scen2_fn << endl;
    scen2.read_data (scen2_fn);

    cout << "Compare" << endl;

    /** @brief Translation map from scenario 1 metabolite indices to scenario 2 metabolite pointers */
    vector<const Metabolite*> xlate_met_12 (scen1.num_metabolites());
    /** @brief Translation map from scenario 2 metabolite indices to scenario 1 metabolite pointers */
    vector<const Metabolite*> xlate_met_21 (scen2.num_metabolites());

    for (auto& met1 : scen1.metabolites()) {
        auto met2 = scen2.find_metabolite (met1->name());
        xlate_met_12[met1->index()] = met2;
        if (met2 != nullptr) {
            if (strcasecmp (
        }
            
        }
    }

    for (auto& react : scen1.reactions()) {
        
    }

    
    return 0;
}
