/*
  String extractor
  - Read agora model dump (0/1 for reaciton in/out of model
  - construct the model
  - find the shortest path 
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
#include "lime/lockfile.h"

#include "mosh/builddate.h" // For _build_date
#include "mosh/scenario.h"
#include "mosh/pathfinder.h"

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

int
main (int argc, const char* argv[]) 
{
    int seed = 0;
    string data_fn = "";
    vector<string> supply_demand_vec;
    string supply_demand_fn = "none";
    string flux_fn = "none";
    string out_fn = "";
    string summary_fn = "";
    string dot_fn = "";
    string debug_str = "";
    string met_name = "";
    bool biomass_mode = false;
    bool explode_react = false;
    bool explode_met = false;

    Config config;
    config.addItem ("seed", seed);
    config.addItem ("data", data_fn);
    config.addItem ("sd_vec", supply_demand_vec);
    config.addItem ("sd_fn", supply_demand_fn);
    config.addItem ("flux_fn", flux_fn);
    

    Opts opts (R"EOF(
    Check path from inputs to biomass reactants
    Loop mode (default): Find loops in reactions
           Biomass mode: Find path to all biomass reactants
    If solution supplied, then only reactions with a non-zero flux in
    that solution are considered
)EOF"
    );

    // Add a switch with a config alternative
    opts.add_opt ("-b", &biomass_mode, "Biomass mode", "biomass");
    opts.add_opt ("-s", &seed, "Rand seed (0 -> current time)", "seed");
    opts.add_opt ("-c", &config, "Config filename");
    opts.add_opt_filename (
        "-sd", &supply_demand_fn, "Supply/demand filename", "sd_fn"
    );
    opts.add_opt_filename (
        "-SD", &supply_demand_vec, "Supply/demand filenames", "sd_vec"
    );
    opts.add_opt_filename ("-f", &flux_fn, "Flux (input sol) filename");
    opts.add_opt_filename ("-o", &out_fn, "Output filename");
    opts.add_opt_filename ("-O", &summary_fn, "Summary filename (appended)");
    opts.add_opt_filename ("-g", &dot_fn, "Solution dot filename");
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

    Rand rand(seed);
    cerr << "Seed " << rand.getSeed() << endl;
    config.addItem ("seed", rand.getSeed());

    Scenario scenario (true, false, 1.0);

    scenario.read_data (data_fn);
    
    if (flux_fn.length() > 0 && flux_fn !="none")
        scenario.read_flux (flux_fn);

    if (supply_demand_fn.length() > 0 && supply_demand_fn != "none") {
        scenario.read_supply_demand (supply_demand_fn);
    }
    else if (supply_demand_vec.size() > 0) {
        for (auto fn : supply_demand_vec) {
            scenario.read_supply_demand (fn);
        }
    }
    else {
        scenario.no_supply_demand();
    }

    std::set<const Metabolite*> dot_mets;
    std::set<const Reaction*> dot_reacts;
    std::set<string> dot_edges;
    
    PathFinder finder (&scenario);

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
    else {
        int num_loops = 0;
        int num_mets = 0;
        if (met_name.length() > 0) {
            auto met = scenario.find_metabolite (met_name);
            if (met == nullptr) {
                opts.usage ("Unknown metabolite " + met_name);
                exit (1);
            }
            while (finder.find_loop (met, num_mets))
                num_loops++;
        }
        else {
            finder.find_loops (num_loops, num_mets);
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
    
    if (out_fn.length() > 0) {
        cout << "Writing " << out_fn << endl;
        ofstream out (out_fn);
        out << "# Produced by " << progname << " on " << todayString() << endl;
        out << "# datafile "  << data_fn << endl;
        out << "# " << summary.str() << endl;
    }
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
