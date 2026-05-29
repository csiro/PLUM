
/*
  LNS style incremental multi-experiment solver
*/

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
#include "lime/linereader.h"
#include "lime/strutil.h"
#include "lime/fileutil.h"
#include "lime/timekeeper.h"
#include "lime/lockfile.h"

#include "mosh/builddate.h" // For _build_date
#include "mosh/scenario.h"
#include "mosh/solution.h"
#include "mosh/lnsmxsolver.h"

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
    vector<string> flavour_name = {"GRB", "LPS", "HIGHS"};
    
    int seed = 0;
    int num_threads = 8;
    int time_limit = 0;
    int max_iters = 1000;
    int which_obj = (int)LnsMxSolver::COST_TAU__ERROR;
    bool by_add = false;
    bool quiet = false;
    bool verbose = false;
    bool make_all_unit = true;
    string flavour_str = flavour_name[0];
    string data_fn = "";
    string react_cost_fn = "none";
    string react_cost_min_fn = "none";
    string sd_list_fn = "";
    string save_best_fn = "none";
    string save_best_sol_fn = "none";
    string progress_fn = "none";
    string out_fn = "";
    string out_sol_fn = "";
    string outmodel_fn = "";
    string summary_fn = "";
    string descr_fn = "";
    string pareto_fn = "";
    string pareto_base_fn = "";
    string debug_str = "";
    
    Config config;
    
    Params params;
    params.init_config (config);

    // Init plummx-only config items
    config.addItem ("flavour", flavour_str);
    config.addItem ("which_obj", which_obj);
    config.addItem ("data", data_fn);
    config.addItem ("react_cost_fn", react_cost_fn);
    config.addItem ("react_cost_min_fn", react_cost_min_fn);
    config.addItem ("sd_lis", sd_list_fn);
    config.addItem ("best_fn", save_best_fn);
    config.addItem ("best_sol_fn", save_best_sol_fn);
    config.addItem ("progress_fn", progress_fn);
    config.addItem ("pareto_fn", pareto_fn);
    config.addItem ("pareto_base_fn", pareto_base_fn);
    
    config.addItem ("iters", max_iters);
    config.addItem ("by_add", by_add);
    config.addItem ("make_all_unit", make_all_unit);
    
    Opts opts (R"EOF(
    Solve the gap-filling problem for microbiome analysis.
)EOF"
    );

    // Add a switch with a config alternative
    opts.add_opt (
        "-s", &params.seed, "Rand seed (0 -> current time)", "seed"
    );
    opts.add_opt (
        "-c", &config, "Config filename"
    );
    opts.add_opt (
        "-A", &by_add, "Solve via an additive alogirhtm", "by_add"
    );
    opts.add_opt (
        "-V", &flavour_str, "Flavour: One of GRB, LPS, HIGHS", "flavour"
    );
    opts.add_opt (
        "-u", &make_all_unit, "Make all reactions unit cost at start",
        "make_all_unit"
    );
    opts.add_opt (
        "-p", &params.num_threads, "Num parallel threads", "threads"
    );
    opts.add_opt (
        "-t", &params.time_limit, "Time limit (secs, 0->none)", "timelimit"
    );
    opts.add_opt (
        "-j", &which_obj, "Which objective (0:error,cost 1:tau+error,cost 2:tau+cost,error 3:cost+tau,error 4:error+cost,tau 5:all 6:tau,cost 7:cost,tau)",
        "which_obj"
    );
    opts.add_opt (
        "-i", &max_iters, "Max LNS iters (0->no limit)", "iters"
    );
    opts.add_opt (
        "-I", &params.num_intervals, "Intervals in target biomass search", "num_intervals"
    );
    opts.add_opt (
        "-F", &params.target_flux,
        "Target flux for rank 0 experiment (0 -> from cts sol)",
        "target_flux"
    );
    opts.add_opt (
        "-T", &params.tau_mult,
        "Multiplier for tau term in objective", "tau_mult"
    );
    opts.add_opt (
        "-CX", &params.max_react_cost,
        "Max react cost - higher is un-selected", "max_react_cost"
    );
    opts.add_opt ("-q", &quiet, "Run quietly", "quiet");
    opts.add_opt ("-v", &verbose, "Verbose");
    opts.add_opt_filename (
        "-rc", &react_cost_fn, "Reaction cost filename", "react_cost_fn"
    );
    opts.add_opt_filename (
        "-RC", &react_cost_min_fn, "Reaction cost filename (min applied)",
        "react_cost_min_fn"
    );
    opts.add_opt_filename (
        "-b", &save_best_fn, "Filename to save intermediate best unselected",
        "best_fn"
    );
    opts.add_opt_filename (
        "-B", &save_best_sol_fn, "Filename to append intermediate best sols",
        "best_sol_fn"
    );
    opts.add_opt_filename (
        "-P", &pareto_fn, "Filename to save pareto front", "pareto_fn"
    );
    opts.add_opt_filename (
        "-PN", &pareto_base_fn, "Base name to save pareto front models (saved as '<fn>-<n>.pld')",
        "pareto_base_fn"
    );
    opts.add_opt_filename (
        "-if", &progress_fn, "Filename to save iteration info",
        "progress_fn"
    );
    opts.add_opt_filename ("-ob", &out_fn, "Output best unselected filename");
    opts.add_opt_filename ("-o", &out_sol_fn, "Output solution filename");
    opts.add_opt_filename ("-M", &outmodel_fn, "Output model filename");
    opts.add_opt_filename ("-O", &summary_fn, "Summary filename (appended)");
    opts.add_opt_filename ("-D", &descr_fn, "Solution description");
    opts.add_opt ("-d", &debug_str, "Debug string");

    // Add a filename args with a config alt
    opts.add_arg ("data.pld", &data_fn, "Input data", "data");
    opts.add_arg (
        "sd.lis", &sd_list_fn, "File of supply/demand filenames", "sd_lis"
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

    TimeKeeper timer;

    assert (data_fn.length() > 0); // Not optional in Opts 

    Rand rand(params.seed);
    DEBUG ('A', "Seed " << rand.getSeed());
    config.addItem ("seed", rand.getSeed());

    params.read_config (config);
    params.finalise();
    cerr << "Params are " << params << endl;
    DEBUG ('A', "Params " << params);

    Scenario scenario;

    scenario.read_data (data_fn);
    
    //scenario.add_dummy_reactions();
    
    LineReader sd_reader (sd_list_fn);
    string line;
    while (sd_reader.getLine (line)) {
        scenario.read_supply_demand (line, &params);
    }
    if (scenario.num_experiments() == 0)
        limeCrash ("No experiments in sd file");

    if (scenario.biomass_react() == nullptr)
        limeCrash ("No biomass reaction in scenario");
    
    if (isFilename (react_cost_fn))
        scenario.read_react_cost (react_cost_fn, &params, mosh::REPLACE);
    if (isFilename (react_cost_min_fn))
        scenario.read_react_cost (react_cost_min_fn, &params, mosh::USEMIN);
    if (params.max_react_cost > 0)
        scenario.unselect_above_cost (params.max_react_cost);
    if (params.unit_cost)
        scenario.set_react_cost (params.gene_ind_cost);
    
    string expstr;
    string sep = "";
    for (auto& exp : scenario.experiments()) {
        expstr += sep + exp->name();
        sep = " ";
    }

    Flavour flavour = GUROBI;
    if (equal_ic (flavour_str, "GRB"))
        flavour = GUROBI;
    else if (equal_ic (flavour_str, "LPS"))
        flavour = LP_SOLVE;
    else if (equal_ic (flavour_str, "HIGHS"))
        flavour = HIGHS;
    else {
        opts.usage ("Unknown falvour: " + flavour_str);
        exit (1);
    }


    LnsMxSolver::WhichObj which_obj_enum = (LnsMxSolver::WhichObj) which_obj;
    LnsMxSolver solver (
        &scenario, &params, flavour, which_obj_enum,
        rand.generateSeed(), max_iters, make_all_unit, progname
    );

    if (isFilename (save_best_fn))
        solver.set_save_best_fn (save_best_fn);
    if (isFilename (save_best_sol_fn)) {
        // Write preamble
        ofstream out (save_best_sol_fn);
        out << "# Produced by " << progname << " on " << todayString() << endl;
        out << "# datafile "  << data_fn << " experiments " << expstr << endl;
        solver.set_save_best_sol_fn (save_best_sol_fn);
    }
    if (isFilename (progress_fn))
        solver.set_progress_fn (progress_fn);
    if (quiet)
        solver.set_quiet (2);
    if (verbose)
        solver.set_quiet (0);
    
    SolutionPtr sol =
        by_add
        ? solver.solve_by_add ()
        : solver.solve();
    auto attribs = solver.best_attr();

    string attrib_str = "FAILED";
    if (attribs != nullptr)
        attrib_str = attribs->to_string();

    stringstream summary;
    summary <<
        "seed " << rand.getSeed() <<
        " data " << data_fn <<
        " experiments " << scenario.num_experiments() <<
        " mets " << scenario.num_metabolites() <<
        " reactions " << scenario.num_reactions() <<
        " numused " << solver.num_used() <<
        " num_nonunit " << solver.num_non_unit() <<
        " " << attrib_str <<
        " iters " << solver.iters() <<
        " numsolves " << solver.inner_iters() <<
        " improves " << solver.num_improves() <<
        " simanneal " << solver.num_simanneal() <<
        " iterfound " << solver.iter_found_best() <<
        " cpu " << solver.elapsed_time_secs() <<
        " wall " << timer.elapsedWallSecs() <<
        " config " << config.show() <<
        " build " << buildId();

    cout << summary.str() << endl;

    stringstream file_header;
    file_header <<
        "# Produced by " << progname << " on " << todayString() << "\n" <<
        "# datafile "  << data_fn << " experiments " << expstr << "\n" <<
        "# " << summary.str();
    
    if (isFilename (out_fn)) {
        cout << "Writing " << out_fn << endl;
        ofstream out (out_fn);
        out << file_header.str() << endl;
        
        if (sol != nullptr)
            solver.write_best_unselected (out);
    }
    if (isFilename (out_sol_fn)) {
        cout << "Writing " << out_sol_fn << endl;
        ofstream out (out_sol_fn);
        out << file_header.str() << endl;
        if (sol != nullptr)
            solver.write_best_sol (out);
    }
    if (isFilename (outmodel_fn)) {
        cout << "Writing " << outmodel_fn << endl;
        ofstream out (outmodel_fn);
        out << file_header.str() << endl;
        solver.write_model (out);
    }
    if (isFilename (pareto_fn)) {
        cout << "Writing " << pareto_fn << endl;
        ofstream out (pareto_fn);
        out << file_header.str() << endl;
        solver.write_pareto (out, pareto_base_fn, file_header.str());
    }
    if (isFilename (summary_fn)) {
        cout << "Writing " << summary_fn << endl;
        LockFile lock (summary_fn);
        ofstream out (summary_fn, fstream::out | fstream::app);
        out << summary.str() << endl;
    }
    if (isFilename (descr_fn)) {
        cout << "Writing " << descr_fn << endl;
        ofstream out (descr_fn);
        out << file_header.str() << endl;
        solver.write_descr (out);
    }
    return 0;
}
