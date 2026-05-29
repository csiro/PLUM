#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>

#include <stdlib.h>

#include "lime/config.h"
#include "lime/opts.h"
#include "lime/rand.h"
#include "lime/dig.h"
#include "lime/debug.h"
#include "lime/error.h"
#include "lime/strutil.h"
#include "lime/linereader.h"
#include "lime/limetok.h"
#include "lime/strutil.h"
#include "lime/lockfile.h"

#include "mosh/builddate.h" // For _build_date
#include "mosh/scenario.h"
#include "mosh/solution.h"
#include "mosh/dummysolver.h"
#include "mosh/mathheur.h" 
#include "mosh/grblpsolver.h" 
#include "mosh/grbintsolver.h" 
#include "mosh/grbcombsolver.h" 

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
    int num_threads = 8;
    int time_limit = 0;
    double biomass_obj_mult = 1.0f;
    double biomass_mult_int = 10.0f;
    int max_biomass_search_iters = 0;
    bool use_abs = false;
    bool use_error = false;
    string data_fn = "";
    string sdfn_fn = "";
    string summary_fn = "";
    string debug_str = "";

    Config config;
    config.addItem ("seed", seed);
    config.addItem ("threads", num_threads);
    config.addItem ("timelimit", time_limit);
    config.addItem ("biomass_mult", biomass_obj_mult);
    config.addItem ("biomass_mult_int", biomass_obj_mult);
    config.addItem ("biomass_iters", max_biomass_search_iters);
    config.addItem ("use_abs", use_abs);
    config.addItem ("use_error", use_error);
    config.addItem ("data", data_fn);
    

    Opts opts (R"EOF(
    Solve the gap-filling problem for microbiome analysis.
)EOF"
    );

    // Add a switch with a config alternative
    opts.add_opt ("-s", &seed, "Rand seed (0 -> current time)", "seed");
    opts.add_opt ("-p", &num_threads, "Num parallel threads", "threads");
    opts.add_opt ("-t", &time_limit, "Time limit (secs, 0->none)", "timelimit");
    opts.add_opt (
        "-a", &use_abs, Opts::TAKES_ARG, "Use absolute objective", "use_abs"
    );
    opts.add_opt (
        "-e", &use_error, Opts::TAKES_ARG, "Use error bounds", "use_error"
    );
    opts.add_opt (
        "-j", &biomass_obj_mult, "Biomass objective multiplier", "biomass_mult"
    );
    opts.add_opt (
        "-ji", &biomass_mult_int, "Biomass objective multiplier", "biomass_mult_int"
    );
    opts.add_opt (
        "-J", &max_biomass_search_iters, "Max iters for positive biomass search", "biomass_iters"
    );
    opts.add_opt ("-c", &config, "Config filename");
    opts.add_opt_filename ("-O", &summary_fn, "Summary filename (appended)");
    opts.add_opt ("-d", &debug_str, "Debug string");

    // Add a filename args with a config alt
    opts.add_arg ("data.pld", &data_fn, "Input data", "data");
    opts.add_arg ("fn.sd", &sdfn_fn, "Supply/demand file of filenames", "sd_fn");
    

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

    vector<string> sd_fn;
    LineReader reader (sdfn_fn);
    string line;
    while (reader.getLine(line)) {
        LimeTok tok (line);
        bool error = false;
        sd_fn.push_back (tok.nextString(error));
    }
    
    Scenario scenario (
        use_abs, use_error, biomass_obj_mult, max_biomass_search_iters,
        biomass_mult_int
    );

    scenario.read_data (data_fn);

    if (flux_fn.length() > 0 && flux_fn !="none")
        scenario.read_flux (flux_fn);

    if (supply_demand_vec.size() > 0) {
        for (auto fn : supply_demand_vec) {
            scenario.read_supply_demand (fn);
        }
    }
    else {
        scenario.no_supply_demand();
    }

    cout << "Setting up solver" << endl;
    GapSolverPtr solver = nullptr;
    try {
        switch (solver_type) {
        case CTS:
            solver =
                make_shared<GrbLPSolver> (
                    &scenario, num_threads, time_limit,
                    rand.generateSeed()
                    
                );
            break;
        case INT:
            solver =
                make_shared<GrbIntSolver> (
                    &scenario, num_threads, time_limit,
                    rand.generateSeed()
                );
            break;
        case DUMMY:
            solver =
                make_shared<DummySolver> (
                    &scenario, rand.generateSeed()
            );
            break;
        case COMB:
            solver =
                make_shared<GrbCombSolver> (
                    &scenario, comb_extra, num_threads, time_limit,
                    rand.generateSeed()
                );
            break;
        case MH:
            solver =
                make_shared<MathHeur> (
                    &scenario, max_mh_iters, num_threads, time_limit,
                    rand.generateSeed()
                );
            break;
        }
    }
    catch(GRBException e) {
        cerr << "Caught exception from Gurobi. Error code = " << e.getErrorCode() << endl;
        cerr << e.getMessage() << endl;
        exit(1);
    }

    cout << "Using solver " << solver_name[solver_type] << endl;

    DigPtr dig = nullptr;
    if (dig_fn.length() > 0) {
        dig = make_shared<Dig> (dig_fn);
        dig->showMessage (string ("Produced by ") + progname);
        dig->showMessage (string ("on ") + todayString());

        Solution dummy  (&scenario);
        dummy.fill (1.0f);
        dummy.draw (dig.get());
        dig->waitAndWipe();
        scenario.draw_reachability (dig.get());
        dig->waitAndWipe();
    }
    
    auto sol = solver->solve();
    
    int reaction_count = 0;
    double sum_flux = 0.0f;
    double obj_value = 0.0f;
    double rel_obj_value = 0.0f;
    double abs_obj_value = 0.0f;
    double biomass_flux = 0.0f;
    double biomass_obj = 0.0f;
    double dummy_flux = 0.0f;
    int sol_depth = 0; 

    if (sol != nullptr) {
        reaction_count = sol->reaction_count();
        sum_flux = sol->sum_flux();
        obj_value = sol->obj_value();
        rel_obj_value = sol->rel_obj_value();
        abs_obj_value = sol->abs_obj_value();
        biomass_flux = sol->biomass_flux();
        biomass_obj = sol->biomass_obj();
        dummy_flux = sol->dummy_flux();
        sol_depth = scenario.sol_depth (sol.get());
    }
                
    int reactions = 0;
    int mets = 0;
    int residuals = 0;
    // Counts the number of reactions, metabolites and metabolites with
    // residuals that *cannot* be reached
    scenario.count_reachability(0, reactions, mets, residuals);
    cout << reactions << " reactions are reachable" << endl;

    stringstream summary;
    summary <<
        "seed " << rand.getSeed() <<
        " data " << data_fn <<
        " mets " << scenario.num_metabolites() <<
        " reactions " << scenario.num_reactions() <<
        " orig-reactions " << scenario.orig_num_reactions() <<
        " numused " << reaction_count <<
        " sumflux " << sum_flux <<
        " obj " << obj_value <<
        " relobj " << rel_obj_value <<
        " absobj " << abs_obj_value <<
        " status " << solver->status_str() << 
        " biomass_flux " << biomass_flux <<
        " biomass_obj " << biomass_obj <<
        " dummy_flux " << dummy_flux <<
        " soldepth " << sol_depth <<
        " depth " << scenario.depth () <<
        " noreach_react " << reactions << 
        " noreach_mets " << mets << 
        " noreach_resid " << residuals <<
        " numselected " << scenario.num_selected() <<
        solver->summary() << 
        " cpu " << solver->elapsed_time_secs() <<
        " config " << config.show() <<
        " build " << buildId();

    cout << summary.str() << endl;
    
    if (out_fn.length() > 0) {
        cout << "Writing " << out_fn << endl;
        ofstream out (out_fn);
        out << "# Produced by " << progname << " on " << todayString() << endl;
        out << "# datafile "  << data_fn << endl;
        out << "# " << summary.str() << endl;

        if (sol != nullptr)
            sol->write_flux (out);
    }
    if (metbal_fn.length() > 0) {
        cout << "Writing " << metbal_fn << endl;
        ofstream out (metbal_fn);
        out << "# Produced by " << progname << " on " << todayString() << endl;
        out << "# datafile "  << data_fn << endl;
        out << "# " << summary.str() << endl;
        
        if (sol != nullptr)
            sol->write_metbal (out);
    }
    if (cs_dot_fn.length() > 0) {
        cout << "Writing " << cs_dot_fn << endl;
        ofstream dot (cs_dot_fn);
        dot << "/* Produced by " << progname << " on " << todayString() <<
            " */"<< endl;
        dot << "/* datafile "  << data_fn << " */" << endl;
        dot << "/* " << summary.str() << " */" << endl;
        
        if (sol != nullptr)
            sol->write_carbon_source_dot (dot);
    }
    if (sd_dot_fn.length() > 0) {
        cout << "Writing " << sd_dot_fn << endl;
        ofstream dot (sd_dot_fn);
        dot << "/* Produced by " << progname << " on " << todayString() <<
            " */"<< endl;
        dot << "/* datafile "  << data_fn << " */" << endl;
        dot << "/* " << summary.str() << " */" << endl;
        
        if (sol != nullptr)
            sol->write_supply_demand_dot (dot);
    }
    if (summary_fn.length() > 0) {
        cout << "Writing " << summary_fn << endl;
        LockFile lock (summary_fn);
        ofstream out (summary_fn, fstream::out | fstream::app);
        out << summary.str() << endl;
    }
    if (dig != nullptr && sol != nullptr) {
        cout << "Writing " << dig_fn << endl;
        
        sol->draw (dig.get());
        dig->waitAndWipe();
        sol->draw_metbal (dig.get());
        dig->waitAndWipe();
        sol->draw_col_val (dig.get());
        dig->waitAndWipe();
        scenario.draw_reachability (dig.get(), sol.get());
    }

    return 0;
}
