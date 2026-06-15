/**
 * @file plumsp.cpp
 * @brief PLUM shortest path analysis tool for metabolic networks
 *
 * This program analyzes metabolic networks to find paths and loops in reactions.
 * It supports two primary modes:
 * - Loop mode (default): Identifies loops in the reaction network
 * - Biomass mode: Finds paths from inputs to all biomass reactants
 *
 * The tool can analyze flux balance analysis (FBA) solutions and consider only
 * reactions with non-zero flux. It supports carbon source graph generation and
 * metabolic gap-filling analysis.
 */
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <set>

#include <stdlib.h>

#include "lime/config.h"
#include "lime/opts.h"
#include "lime/rand.h"
#include "lime/dig.h"
#include "lime/debug.h"
#include "lime/error.h"
#include "lime/strutil.h"
#include "lime/fileutil.h"
#include "lime/linereader.h"
#include "lime/lockfile.h"

#include "mosh/builddate.h" // For _build_date
#include "mosh/scenario.h"
#include "mosh/params.h"
#include "mosh/pathfinder.h"

using namespace std;
using namespace lime;
using namespace mosh;

// Define the OS variable
/**
 * @brief Determines the operating system at compile time
 *
 * Uses preprocessor macros to identify the target operating system
 * during compilation.
 *
 * @return Constant string identifying the OS ("Cygwin", "Windows", "Mac OS X", "Linux", "Unix", or "Unknown OS")
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
 * @brief Main entry point for the PLUM shortest path analysis tool
 *
 * Parses command-line arguments, initializes the metabolic network scenario,
 * and performs path/loop analysis based on the selected mode:
 * - Loop mode: Identifies reaction loops in the metabolic network
 * - Biomass mode: Analyzes paths to biomass reactants
 * - Carbon depth mode: Generates carbon source graphs
 *
 * The program can optionally read flux solutions, reaction costs, and
 * supply/demand constraints to refine the analysis. Results can be written
 * to output files in various formats including DOT graphs for visualization.
 *
 * @param argc Number of command-line arguments
 * @param argv Array of command-line argument strings
 * @return Exit status (0 for success, 1 for failure)
 */
int
main (int argc, const char* argv[])
{
    int seed = 0; /**< Random number generator seed (0 = use current time) */
    string data_fn = ""; /**< Input data filename (metabolic network file) */
    string supply_demand_lis = ""; /**< File containing list of supply/demand filenames */ 
    string supply_demand_fn = "none"; /**< Supply/demand constraints filename */
    string flux_fn = "none"; /**< Flux solution input filename (FBA solution) */
    string react_cost_fn = "none"; /**< Reaction cost filename for gap-filling optimization */
    string react_cost_min_fn = "none"; /**< Minimum reaction cost filename */
    string out_fn = ""; /**< Output filename for analysis results */
    string summary_fn = ""; /**< Summary filename (results appended) */
    string dot_fn = ""; /**< DOT format graph output filename for visualization */
    string debug_str = ""; /**< Debug mode filter string */
    string met_name = ""; /**< Metabolite name for single-metabolite analysis */
    int carbon_depth = 0; /**< Depth for carbon source graph traversal */
    bool biomass_mode = false; /**< Enable biomass mode (find paths to biomass reactants) */
    bool explode_react = false; /**< Expand reactions in DOT graph visualization */
    bool explode_met = false; /**< Expand metabolites in DOT graph visualization */
    double min_flux = 0.0f; /**< Minimum flux threshold for reporting reactions */

    Config config; /**< Configuration object for parameter management */

    Params params; /**< Parameters object for metabolic network analysis */
    params.init_config (config);
    
    config.addItem ("seed", seed);
    config.addItem ("data", data_fn);
    config.addItem ("sd_lis", supply_demand_lis);
    config.addItem ("sd_fn", supply_demand_fn);
    config.addItem ("flux_fn", flux_fn);
    config.addItem ("react_cost_fn", react_cost_fn);
    config.addItem ("react_cost_min_fn", react_cost_min_fn);

    Opts opts (R"EOF(
    Check path from inputs to biomass reactants
    Loop mode (default): Find loops in reactions
           Biomass mode: Find path to all biomass reactants
    If solution supplied, then only reactions with a non-zero flux in
    that solution are considered
)EOF"
    ); /**< Command-line options parser with usage description */

    // Add a switch with a config alternative
    opts.add_opt ("-b", &biomass_mode, "Use biomass mode", "biomass");
    opts.add_opt ("-C", &carbon_depth, "Draw carbon source graph - arg gives depth", "carbon_depth");
    opts.add_opt ("-s", &seed, "Rand seed (0 -> current time)", "seed");
    opts.add_opt ("-c", &config, "Config filename");
    opts.add_opt_filename (
        "-sd", &supply_demand_fn, "Supply/demand filename", "sd_fn"
    );
    opts.add_opt_filename (
        "-SD", &supply_demand_lis, "File of supply/demand filenames", "sd_lis"
    );
    opts.add_opt_filename (
        "-rc", &react_cost_fn, "Reaction cost filename", "react_cost_fn"
    );
    opts.add_opt (
        "-CX", &params.max_react_cost,
        "Max react cost - higher is un-selected", "max_react_cost"
    );
    opts.add_opt_filename ("-f", &flux_fn, "Flux (input sol) filename");
    opts.add_opt_filename ("-o", &out_fn, "Output filename");
    opts.add_opt_filename ("-O", &summary_fn, "Summary filename (appended)");
    opts.add_opt_filename ("-g", &dot_fn, "Dot filename");
    opts.add_opt ("-F", &min_flux, "Min flux to report");
    opts.add_opt ("-x", &explode_react, "Explode reactions");
    opts.add_opt ("-X", &explode_met, "Explode metabolites");
    opts.add_opt ("-d", &debug_str, "Debug string");

    // Add a filename args with a config alt
    opts.add_arg ("data.pld", &data_fn, "Input data", "data");
    
    opts.add_optional_arg (
        "met", &met_name, "Metabolite name (for single-met run)"
    );

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

    if (!opts.process (argc, argv, &config))
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

    Rand rand(seed); /**< Random number generator initialized with seed */
    cerr << "Seed " << rand.getSeed() << endl;
    config.addItem ("seed", rand.getSeed());

    Scenario scenario; /**< Metabolic network scenario containing reactions and metabolites */

    scenario.read_data (data_fn);
    
    if (flux_fn.length() > 0 && flux_fn !="none")
        scenario.read_flux (flux_fn);

    if (isFilename (react_cost_fn))
        scenario.read_react_cost (react_cost_fn, &params, mosh::REPLACE);
    if (isFilename (react_cost_min_fn))
        scenario.read_react_cost (react_cost_min_fn, &params, mosh::USEMIN);
    if (params.max_react_cost > 0)
        scenario.unselect_above_cost (params.max_react_cost);

    if (supply_demand_fn.length() > 0 && supply_demand_fn != "none") {
        scenario.read_supply_demand (supply_demand_fn, &params);
    }
    else if (supply_demand_lis.size() > 0) {
        LineReader reader (supply_demand_lis);
        string fn;
        while (reader.getLine (fn)) {
            scenario.read_supply_demand (fn, &params);
        }
    }
    else {
        scenario.no_supply_demand();
    }

    std::set<const Metabolite*> dot_mets; /**< Set of metabolites for DOT graph output */
    std::set<const Reaction*> dot_reacts; /**< Set of reactions for DOT graph output */
    std::set<string> dot_edges; /**< Set of edges for DOT graph output */

    PathFinder finder (&scenario, &params); /**< PathFinder object for analyzing paths and loops in the metabolic network */

    if (out_fn.length() > 0) {
        ofstream_ptr out = make_shared<ofstream> (out_fn);
        *out << "# Produced on " << todayString() << " by " << progname << endl;
        *out << "# datafile "  << data_fn << endl;
        finder.set_outfile (out);
    }
    
    stringstream mode_summary;
    
    if (biomass_mode) {
        int count = 0;
        int count_ok = 0;
        finder.check_biomass (count, count_ok);
        
        mode_summary <<
            " biomass_mets " << count <<
            " count_ok " << count_ok <<
            " ok_pc " << (int) round ((100.0f * count_ok) / count);
    }
    else if (carbon_depth > 0) {
        if (!isFilename (dot_fn))
            opts.usage ("Dot filename required for carbon source graph");
        if (met_name.length() == 0) 
            opts.usage ("Met name required for carbon source graph");
        auto met = scenario.find_metabolite (met_name);
        if (met == nullptr) 
            opts.usage ("Unknown metabolite " + met_name);
        ofstream dot (dot_fn);
        dot << "// Produced on " + todayString() + " by " + progname << endl;
        
        finder.carbon_source_graph (met, carbon_depth, dot);
    }
    else {
        int num_loops = 0;
        int num_mets = 0;
        if (met_name.length() > 0) {
            auto met = scenario.find_metabolite (met_name);
            if (met == nullptr) 
                opts.usage ("Unknown metabolite " + met_name);
            while (finder.find_loop (met, min_flux, num_mets))
                num_loops++;
        }
        else {
            finder.find_loops (min_flux, num_loops, num_mets);
        }
        mode_summary <<
            " num_loops " << num_loops <<
            " num_mets_involved " << num_mets;
    }

    stringstream summary;
    summary <<
        "seed " << rand.getSeed() <<
        " data " << data_fn <<
        " mets " << scenario.num_metabolites() <<
        " reactions " << scenario.num_reactions() <<
        " edges " << finder.num_edges() <<
        mode_summary.str() <<
        " config " << config.show() <<
        " build " << buildId();

    cout << summary.str() << endl;
    
    if (summary_fn.length() > 0) {
        cout << "Writing " << summary_fn << endl;
        LockFile lock (summary_fn);
        ofstream out (summary_fn, fstream::out | fstream::app);
        out << summary.str() << endl;
    }
    if (dot_fn.length() > 0) {
        cout << "Writing " << dot_fn << endl;
        ofstream dot (dot_fn);
        finder.write_dot (dot, explode_react, explode_met);
        cout << "View with" << endl;
        cout << "xdot " << dot_fn << " &" << endl;
    }
    
    return 0;
}
