/**
 * @file plum.cpp
 * @brief Main entry point for the PLUM metabolic gap-filling solver
 *
 * This file implements the PLUM (Probabilistic modeling for metabolomics) application,
 * which solves metabolic gap-filling problems for microbiome analysis using various
 * optimization approaches including continuous (CTS), integer (INT), and incremental
 * (INCR) solvers with flux balance analysis (FBA) techniques.
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
#include "lime/strutil.h"
#include "lime/fileutil.h"
#include "lime/linereader.h"
#include "lime/lockfile.h"
#include "lime/numutil.h"

#include "mosh/builddate.h" // For _build_date
#include "mosh/constants.h"
#include "mosh/scenario.h"
#include "mosh/params.h"
#include "mosh/solution.h"
#include "mosh/dummysolver.h"
#include "mosh/lpsolver.h" 
#ifdef PLUM_GUROBI
#include "mosh/grblpsolverimp.h"
#endif
#ifdef PLUM_LPSOLVER
#include "mosh/lpslpsolverimp.h" 
#endif
#ifdef PLUM_HIGHS
#include "mosh/highslpsolverimp.h" 
#endif
#include "mosh/intsolver.h" 
#include "mosh/incrsolver.h" 

using namespace std;
using namespace lime;
using namespace mosh;

// Define the OS variable
/**
 * @brief Detects and returns the operating system name
 * @return C-string representing the current operating system (e.g., "Linux", "Mac OS X", "Windows")
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
 * @brief Specifies which set of reactions to include in output model
 * @enum WhichReactions
 */
/** @var WhichReactions::USED
 * Only reactions with non-zero flux in the solution
 */
/** @var WhichReactions::GENE_IND
 * Gene-indicated reactions including those with minimum cost
 */
/** @var WhichReactions::SELECTED
 * All selected reactions regardless of flux
 */
enum WhichReactions {USED, GENE_IND, SELECTED};
/**< String names corresponding to WhichReactions enum values */
vector<const char*> whichName = {"USED", "GENE_IND", "SELECTED"};

/**
 * @brief Writes a metabolic model to file with optional reactions included
 * @param fn Output filename for the model
 * @param which Specifies which reactions to include (USED, GENE_IND, or SELECTED)
 * @param react_cost_one If true, set all reaction costs to 1 in output
 * @param scenario Pointer to the metabolic scenario containing network data
 * @param params Pointer to parameter configuration
 * @param sol Shared pointer to the solution containing flux values
 * @param progname Program name string for file header
 * @param summary Summary string with run statistics for file header
 */
void
write_model (
    string fn, WhichReactions which, bool react_cost_one,
    const Scenario* scenario, const Params* params, SolutionPtr sol,
    string progname, string summary
)
{
    cout << "Writing " << fn << endl;
    ofstream out (fn);
    out << "# Produced by " << progname << " on " << todayString() << endl;
    out << "# Which reactions " << whichName[which] <<
        " set all react costs to 1? " << react_cost_one << endl;
    out << "# " << summary << endl;
        
    if (sol != nullptr) {
        sol->write_model (out, react_cost_one);

        if (which == GENE_IND) {
            // Add all reacts with cost min-cost that are not in the model
            out << "# MinCost reactions" << endl;
            for (auto& react_ptr : scenario->reactions()) {
                auto react = react_ptr.get();
                if (
                    params->is_gene_indicated(react) &&
                    limeIsZero (sol->flux (react))
                ) {
                    react->write_to (out);
                }
            }
        }
        else if (which == SELECTED) {
            // Add all selected reacts with no flux
            out << "# MinCost reactions" << endl;
            for (auto& react_ptr : scenario->reactions()) {
                auto react = react_ptr.get();
                if (
                    react->is_selected() &&
                    limeIsZero (sol->flux (react))
                ) {
                    react->write_to (out);
                }
            }
        }
    }
}

/**
 * @brief Main entry point for PLUM gap-filling solver
 *
 * Parses command-line arguments, reads metabolic network data, configures and runs
 * the selected solver (CTS, INT, INCR, etc.), and outputs results including flux
 * distributions, models, and metabolite balances.
 *
 * @param argc Argument count
 * @param argv Argument vector containing command-line parameters
 * @return Exit status (0 for success, 1 for error)
 */
int
main (int argc, const char* argv[])
{
    /**
     * @brief Available solver types for metabolic gap-filling optimization
     * @enum SolverType
     */
    /** @var SolverType::CTS
     * Continuous solver (linear programming)
     */
    /** @var SolverType::CTS2
     * Two-stage continuous solver with target flux
     */
    /** @var SolverType::INT
     * Integer solver (mixed-integer programming)
     */
    /** @var SolverType::INT2
     * Integer solver variant 2
     */
    /** @var SolverType::INT3
     * Integer solver variant 3
     */
    /** @var SolverType::INCR
     * Incremental solver that progressively adds reactions
     */
    /** @var SolverType::INCR_INT
     * Incremental solver with integer programming
     */
    /** @var SolverType::COMB
     * Combinatorial solver (not currently implemented)
     */
    /** @var SolverType::MH
     * Math-heuristic solver (not currently implemented)
     */
    /** @var SolverType::DUMMY
     * Dummy solver for testing purposes
     */
    enum SolverType {
        CTS, CTS2, INT, INT2, INT3, INCR, INCR_INT, COMB, MH, DUMMY
    };
    vector<string> solver_name =
        {
            "CTS", "CTS2", "INT", "INT2", "INT3", "INCR", "INCR_INT", "COMB",
            "MH", "DUMMY"
        };
    vector<string> flavour_name = {"GRB", "LPS", "HIGHS"};
    
    bool quiet = false;
    string solver_str = solver_name[0];
    string flavour_str = flavour_name[0];
    
    string data_fn = "";
    string flux_fn = "none";
    string react_cost_fn = "none";
    string react_cost_min_fn = "none";
    string supply_demand_lis = "none";
    string supply_demand_fn = "none";
    string c_source_fn = "none";
    string cycle_met_fn = "none";
    string base_flux_fn = "none";
    string out_fn = "";
    string cts_out_fn = "";
    string outmodel_fn = "";
    string outmodel1_fn = "";
    string outmodel_plus_fn = "";
    string outmodel_all_fn = "";
    string cs_dot_fn = "";
    string sd_dot_fn = "";
    string met_bal_fn = "";
    string met_use_fn = "";
    string met_dual_fn = "";
    string selected_react_fn = "";
    string bestsol_fn = "";
    string summary_fn = "";
    string dig_fn = "";
    string debug_str = "";
    double min_cost = 1.0;

    Config config;
    DEBUG ('C',"Config " << config);

    Params params;
    params.init_config (config);

    // Init plum-only config items
    config.addItem ("solver", solver_str);
    config.addItem ("flavour", flavour_str);
    
    config.addItem ("data", data_fn);
    config.addItem ("flux_fn", flux_fn);
    config.addItem ("react_cost_fn", react_cost_fn);
    config.addItem ("react_cost_min_fn", react_cost_min_fn);
    config.addItem ("sd_lis", supply_demand_lis);
    config.addItem ("sd_fn", supply_demand_fn);
    config.addItem ("c_source_fn", c_source_fn);
    config.addItem ("cycle_met_fn", cycle_met_fn);
    config.addItem ("base_flux_fn", base_flux_fn);

    Opts opts (R"EOF(
    Solve the gap-filling problem for microbiome analysis.
)EOF"
    );

    // Add a switch with a config alternative
    opts.add_opt (
        "-s", &params.seed, "Rand seed (0 -> current time)", "seed"
    );
    opts.add_opt (
        "-p", &params.num_threads, "Num parallel threads", "threads"
    );
    opts.add_opt (
        "-t", &params.time_limit, "Time limit (secs, 0->none)", "timelimit"
    );
    opts.add_opt (
        "-v", &solver_str,
        "Solver: One of CTS, CTS2, INT, INT2, INT3, INCR, INCR_INT, COMB, or DUMMY",
        "solver"
    );
    opts.add_opt (
        "-V", &flavour_str, "Flavour: One of GRB, LPS, HIGHS", "flavour"
    );
    opts.add_opt (
        "-a", &params.use_abs_obj, Opts::TAKES_ARG,
        "Use absolute objective", "use_abs"
    );
    opts.add_opt (
        "-e", &params.use_error, Opts::TAKES_ARG,
        "Use error bounds", "use_error"
    );
    opts.add_opt (
        "-D", &params.use_dummy, "Use dummy reactions", "use_dummy"
    );
    opts.add_opt (
        "-DRD", &params.use_dummy_biomass_react_dm,
        "Use dummy DM reactions for biomass Reactants", "use_dummy_biomass_react_dm"
    );
    opts.add_opt (
        "-DRX", &params.use_dummy_biomass_react_ex,
        "Use dummy EX reactions for biomass Reactants", "use_dummy_biomass_react_ex"
    );
    opts.add_opt (
        "-DPD", &params.use_dummy_biomass_prod_dm,
        "Use dummy DM reactions for biomass Products", "use_dummy_biomass_prod_dm"
    );
    opts.add_opt (
        "-DK", &params.preserve_dummies,
        "Keep dummy reactions (even expensive ones)", "preserve_dummies"
    );
    opts.add_opt (
        "-B", &params.biomass_ub, "Biomass flux upper bound (0 -> use default)",
        "biomass_ub"
    );
    opts.add_opt (
        "-ra", &params.detect_runaway,
        "Turn off runaway detection", "detect_runaway"
    );
    opts.add_opt (
        "-q", &quiet, "Quieter run"
    );
    opts.add_opt (
        "-w", &params.biomass_wiggle_pc,
        "Wiggle room for biomass metabolites (as a %)",
        "biomass_wiggle_pc"
    );
    opts.add_opt (
        "-C", &min_cost,
        "Min acceptable cost in Incremental solver", "min_cost"
    );
    opts.add_opt (
        "-CX", &params.max_react_cost,
        "Max react cost - higher is un-selected", "max_react_cost"
    );
    opts.add_opt (
        "-1", &params.unit_cost,
        "Unit cost - all reactions (except biomass) cost 1", "unit_cost"
    );
    opts.add_opt (
        "-j", &params.init_biomass_obj_mult,
        "Biomass objective multiplier", "init_biomass_mult"
    );
    opts.add_opt (
        "-ji", &params.biomass_mult_int,
        "2nd objective multiplier for int solver", "biomass_mult_int"
    );
    opts.add_opt (
        "-J", &params.max_biomass_search_iters,
        "Max iters for positive biomass search", "max_biomass_iters"
    );
    opts.add_opt (
        "-h", &params.max_mh_iters, "Maximum number of math-heuristic iters",
        "max_mh_iters"
    );
    opts.add_opt (
        "-c", &config, "Config filename"
    );
    opts.add_opt_filename (
        "-f", &flux_fn, "Known flux filename", "flux_fn"
    );
    opts.add_opt (
        "-F", &params.target_flux,
        "Target flux for rank 0 experiment (0 -> from cts sol)", "target_flux"
    );
    opts.add_opt_filename (
        "-rc", &react_cost_fn, "Reaction cost filename", "react_cost_fn"
    );
    opts.add_opt_filename (
        "-RC", &react_cost_min_fn, "Reaction cost filename (min applied)",
        "react_cost_min_fn"
    );
    opts.add_opt_filename (
        "-sd", &supply_demand_fn, "Supply/demand filename", "sd_fn"
    );
    opts.add_opt_filename (
        "-SD", &supply_demand_lis, "File of supply/demand filenames", "sd_lis"
    );
    opts.add_opt_filename (
        "-cs", &c_source_fn, "Carbon source metabolites filename", "c_source_fn"
    );
    opts.add_opt_filename (
        "-cm", &cycle_met_fn, "Cycle-met filename", "cycle_met_fn"
    );
    opts.add_opt_filename (
        "-bf", &base_flux_fn, "Base flux from cts sol filename ",
        "base_flux_fn"
    );
    opts.add_opt_filename (
        "-o", &out_fn, "Output filename"
    );
    opts.add_opt_filename (
        "-oc", &cts_out_fn, "Output filename (CTS sols for INT solver)"
    );
    opts.add_opt_filename (
        "-M", &outmodel_fn, "Output model filename"
    );
    opts.add_opt_filename (
        "-M1", &outmodel1_fn, "Output model filename (react cost 1)"
    );
    opts.add_opt_filename (
        "-M+", &outmodel_plus_fn, "Output model filename incl. react cost 1"
    );
    opts.add_opt_filename (
        "-MA", &outmodel_all_fn, "Output model filename - all selected reactions"
    );
    opts.add_opt_filename (
        "-mb", &met_bal_fn, "Metabolite balance filename"
    );
    opts.add_opt_filename (
        "-mu", &met_use_fn, "Metabolite use filename"
    );
    opts.add_opt_filename (
        "-md", &met_dual_fn, "Metabolite dual vals filename"
    );
    opts.add_opt_filename (
        "-rs", &selected_react_fn, "Selected reaction filename"
    );
    opts.add_opt_filename (
        "-gv", &cs_dot_fn, "Carbon-source graph filename"
    );
    opts.add_opt_filename (
        "-gvsd", &sd_dot_fn, "Supply/demand graph filename"
    );
    opts.add_opt_filename (
        "-O", &summary_fn, "Summary filename (appended)"
    );
    opts.add_opt_filename (
        "-b", &bestsol_fn, "Best sol filename ('_<n>.sol' appended)"
    );
    opts.add_opt_filename (
        "-g", &dig_fn, "Dig filename"
    );
    opts.add_opt ("-d", &debug_str, "Debug string");

    // Add a filename args with a config alt
    opts.add_arg ("data.pld", &data_fn, "Input data", "data");
    

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
    cout << "Debug " << debug_str << endl;
    Debug::setFilename ("debug.out");
    Debug::setKey (debug_str);
    Debug::debugFile() << progname << endl;
    Debug::debugFile() << "Timestamp " << lime::todayString() << endl;
    Debug::debugFile() << "DebugStr: " << debug_str << endl;
#endif

    TimeKeeper timer;

    cout << "Config " << config << endl;
    DEBUG ('A', "Config is " << config); 
    
    params.read_config (config);
    
    assert (data_fn.length() > 0); // Not optional in Opts

    Rand rand(params.seed);
    cerr << "Seed " << rand.getSeed() << endl;
    config.addItem ("seed", rand.getSeed());

    SolverType solver_type = CTS;
    if (equal_ic (solver_str, "CTS"))
        solver_type = CTS;
    else if (equal_ic (solver_str, "CTS2"))
        solver_type = CTS2;
    else if (equal_ic (solver_str, "INT"))
        solver_type = INT;
    else if (equal_ic (solver_str, "INT2"))
        solver_type = INT2;
    else if (equal_ic (solver_str, "INT3"))
        solver_type = INT3;
    else if (equal_ic (solver_str, "INCR"))
        solver_type = INCR;
    else if (equal_ic (solver_str, "INCR_INT"))
        solver_type = INCR_INT;
    else if (equal_ic (solver_str, "COMB"))
        solver_type = COMB;
    else if (equal_ic (solver_str, "MH"))
        solver_type = MH;
    else if (equal_ic (solver_str, "DUMMY"))
        solver_type = DUMMY;
    else {
        opts.usage ("Unknown solver: " + solver_str);
        exit (1);
    }
    if (
        solver_type == INT || solver_type == INT2 || solver_type == INT3 ||
        solver_type == COMB
    ) {
        params.use_abs_obj = true;
        config.addItem ("use_abs", params.use_abs_obj);
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
    
    params.finalise();
    DEBUG ('A', "Params are " << params); 

    
    Scenario scenario;

    scenario.read_data (data_fn);
    if (scenario.biomass_react() == nullptr)
        limeCrash ("No biomass reaction in scenario");

    if (isFilename (flux_fn))
        scenario.read_flux (flux_fn);

    if (isFilename (react_cost_fn))
        scenario.read_react_cost (react_cost_fn, &params, mosh::REPLACE);
    if (isFilename (react_cost_min_fn))
        scenario.read_react_cost (react_cost_min_fn, &params, mosh::USEMIN);
        
    if (isFilename (cycle_met_fn))
        scenario.read_cycle_mets (cycle_met_fn);
    if (isFilename (c_source_fn))
        scenario.read_c_sources (c_source_fn);

    if (isFilename (supply_demand_fn))
        scenario.read_supply_demand (supply_demand_fn, &params);
    else if (isFilename (supply_demand_lis)) {
        LineReader reader (supply_demand_lis);
        string fn;
        while (reader.getLine (fn)) {
            scenario.read_supply_demand (fn, &params);
        }
    }
    else {
        scenario.no_supply_demand();
    }

    if (isFilename (base_flux_fn))
        scenario.read_base_flux (base_flux_fn);
    
    scenario.finalise (&params);

    ofstream selected_react_file;
    if (isFilename (selected_react_fn)) {
        cout << "Opening " << selected_react_fn << endl;
        selected_react_file.open (selected_react_fn);
        selected_react_file <<
            "# Produced by " << progname << " on " << todayString() << endl;
        selected_react_file << "# datafile "  << data_fn << endl;
    }

    cout << "Setting up solver" << endl;
    int imp_seed = rand.uniform0n_1 (1000000);
    int imp_seed2 = rand.uniform0n_1 (1000000);
    GrbEnvPtr grbenv = nullptr;
    LPSolverImpPtr imp = nullptr;
    LPSolverImpPtr imp2 = nullptr;
    bool need_imp2 =
        solver_type == INT ||
        solver_type == INT2 ||
        solver_type == INT3 ||
        solver_type == INCR_INT;

    switch (flavour) {
    case GUROBI:
        try {
            grbenv = make_shared<GRBEnv>();
        }
        catch (GRBException e) {
            limeCrash (
                "Caught Gurobi exception " << e.getErrorCode() << ": " <<
                e.getMessage()
            );
        }
        imp = make_shared<GrbLPSolverImp> (
            &scenario, &params, imp_seed, grbenv
        );
        if (need_imp2) {
            imp2 =
                make_shared<GrbLPSolverImp> (
                    &scenario, &params, imp_seed2, grbenv
                );
        }
        break;

    case LP_SOLVE:
#ifdef PLUM_LPSOLVER
        imp = make_shared<LpsLPSolverImp> (&scenario, &params);
        if (need_imp2) {
            imp2 = make_shared<LpsLPSolverImp> (&scenario, &params);
        }
#else
        limeCrash ("LPSOLVER not compiled in");
#endif
        break;
        
    case HIGHS:
#ifdef PLUM_HIGHS
        imp = make_shared<HighsLPSolverImp> (&scenario, &params, imp_seed);
        if (need_imp2) {
            imp2 = make_shared<HighsLPSolverImp> (&scenario, &params, imp_seed);
        }
#else
        limeCrash ("HIGHS solver not compiled in");
#endif
        break;
    }
    assert (imp != nullptr);
    assert (!need_imp2 || imp2 != nullptr);

    MultiSolPtr cts_sol = nullptr;
    if (
        solver_type == INT ||
        solver_type == INT2 ||
        solver_type == INT3
    ) {
        // Need continuous solutions. Generate them here.
        DEBUG ('A', "Solve experiments using CTS");
        cts_sol = make_shared<MultiSol> (&scenario, &params);
        LPSolver cts_solver (
            &scenario, &params, rand.generateSeed(), 1, imp2
        );
        for (
            size_t exp_idx = 0; exp_idx < scenario.num_experiments(); exp_idx++
        ) {
            const Experiment* exp = scenario.experiment(exp_idx);
            DEBUG ('A', "  Solve for exp " << exp->name());
            cout << "Find CTS sol for " << exp->name() << endl;
            cts_solver.set_which_experiment (exp_idx);
            SolutionPtr exp_sol = cts_solver.solve();
            if (exp_sol == nullptr) {
                DEBUG ('A', "Continuous solve failed");
                limeCrash ("Continuous solve failed for " + exp->name());
            }
            else {
                DEBUG ('A', "  biomass flux is " << exp_sol->biomass_flux());
                cerr << "  Biomass flux for " << exp->name() <<
                    " is " << exp_sol->biomass_flux() << endl;
            }
            cts_sol->set_sol (exp_idx, exp_sol);
        }
        // Release the imp2 memory
        imp2 = nullptr;
    }
    
    GapSolverPtr solver = nullptr;
    switch (solver_type) {
    case CTS:
        solver =
            make_shared<LPSolver> (
                &scenario, &params, rand.generateSeed(), 1, imp
            );
        break;
    case CTS2: {
        auto solver2 =
            make_shared<LPSolver> (
                &scenario, &params, rand.generateSeed(), 2, imp
            );
        solver2->set_target_flux (params.target_flux);
        params.max_biomass_search_iters = 0;
        solver = solver2;
        break;
    }
    case INT: 
    case INT2: 
    case INT3: 
        solver =
            make_shared<IntSolver> (
                &scenario, &params, rand.generateSeed(),
                (solver_type - INT) + 1, imp, cts_sol
            );
    
        break;
    case INCR:
        solver =
            make_shared<IncrSolver> (
                &scenario, &params, rand.generateSeed(),
                imp, nullptr,
                min_cost, selected_react_file
            );
        break;
    case INCR_INT:
        solver =
            make_shared<IncrSolver> (
                &scenario, &params, rand.generateSeed(),
                imp, imp2, 
                min_cost, selected_react_file
            );
        break;
    case DUMMY:
        solver =
            make_shared<DummySolver> (
                &scenario, &params, rand.generateSeed()
            );
        break;
    case COMB:
        /*
          solver =
          make_shared<GrbCombSolver> (
          &scenario, &params, grbenv, 
          rand.generateSeed()
          );
        */
        break;
    case MH:
        /*
          solver =
          make_shared<MathHeur> (
          &scenario, &params, grbenv, rand.generateSeed()
          );
        */
        break;
    }

    cout << "Using solver " << solver_name[solver_type] << endl;

    if (isFilename (bestsol_fn)) 
        solver->set_best_sol_fn (bestsol_fn);

    DigPtr dig = nullptr;
    if (isFilename (dig_fn)) {
        dig = make_shared<Dig> (dig_fn);
        dig->showMessage (string ("Produced by ") + progname);
        dig->showMessage (string ("on ") + todayString());

        Solution dummy  (&scenario, &params);
        dummy.fill (1.0f);
        dummy.draw (dig.get());
        dig->waitAndWipe();
        scenario.draw_reachability (dig.get());
        dig->waitAndWipe();
    }
    solver->set_quiet(quiet);
    
    auto sol = solver->solve();
    
    int reaction_count = 0;
    double sum_flux = 0.0f;
    double rel_obj_value = 0.0f;
    double abs_obj_value = 0.0f;
    double biomass_flux = 0.0f;
    double biomass_obj = 0.0f;
    double dummy_flux = 0.0f;
    double dummy_obj = 0.0f;
    double max_cost = 0.0f;
    int num_dummy = 0.0f;
    int sol_depth = 0; 

    if (sol != nullptr) {
        reaction_count = sol->reaction_count();
        sum_flux = sol->sum_flux();
        rel_obj_value = sol->rel_obj_value();
        abs_obj_value = sol->abs_obj_value();
        biomass_flux = sol->biomass_flux();
        dummy_flux = sol->dummy_flux();
        dummy_obj = sol->dummy_obj_val();
        max_cost = sol->max_cost();
        num_dummy = sol->num_dummy();
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
        " relobj " << rel_obj_value <<
        " absobj " << abs_obj_value <<
        " dummyobj " << dummy_obj <<
        " max_used_cost " << max_cost <<
        " status " << solver->status_str() << 
        " biomass_flux " << biomass_flux <<
        " dummy_flux " << dummy_flux <<
        " num_dummy " << num_dummy <<
        " soldepth " << sol_depth <<
        " depth " << scenario.depth () <<
        " noreach_react " << reactions << 
        " noreach_mets " << mets << 
        " noreach_resid " << residuals <<
        " numselected " << scenario.num_selected() <<
        solver->summary() << 
        " cpu " << solver->elapsed_time_secs() <<
        " wall " << timer.elapsedWallSecs() <<
        " config " << config.show() <<
        " build " << buildId();

    cout << summary.str() << endl;
    
    if (isFilename (out_fn)) {
        cout << "Writing " << out_fn << endl;
        ofstream out (out_fn);
        out << "# Produced by " << progname << " on " << todayString() << endl;
        out << "# datafile "  << data_fn << endl;
        out << "# " << summary.str() << endl;

        if (sol != nullptr)
            sol->write_flux (out);
    }
    if (isFilename (cts_out_fn)) {
        cout << "Writing " << cts_out_fn << endl;
        ofstream out (cts_out_fn);
        out << "# Produced by " << progname << " on " << todayString() << endl;
        out << "# datafile "  << data_fn << endl;
        out << "# " << summary.str() << endl;

        if (cts_sol != nullptr)
            cts_sol->write_flux (out);
    }
    if (isFilename (outmodel_fn)) {
        write_model (
            outmodel_fn, USED, false, &scenario, &params, sol,
            progname, summary.str()
        );
    }
    if (isFilename (outmodel1_fn)) {
        write_model (
            outmodel1_fn, USED, true, &scenario, &params, sol,
            progname, summary.str()
        );
    }
    if (isFilename (outmodel_plus_fn)) {
        write_model (
            outmodel_plus_fn, GENE_IND, true, &scenario, &params, sol,
            progname, summary.str()
        );
    }
    if (isFilename (outmodel_all_fn)) {
        write_model (
            outmodel_all_fn, SELECTED, false, &scenario, &params, sol,
            progname, summary.str()
        );
    }
    if (isFilename (met_bal_fn)) {
        cout << "Writing " << met_bal_fn << endl;
        ofstream out (met_bal_fn);
        out << "# Produced by " << progname << " on " << todayString() << endl;
        out << "# datafile "  << data_fn << endl;
        out << "# " << summary.str() << endl;
        
        if (sol != nullptr)
            sol->write_metbal (out);
    }
    if (isFilename (met_use_fn)) {
        cout << "Writing " << met_use_fn << endl;
        ofstream out (met_use_fn);
        out << "# Produced by " << progname << " on " << todayString() << endl;
        out << "# datafile "  << data_fn << endl;
        out << "# " << summary.str() << endl;
        
        if (sol != nullptr)
            sol->write_met_use (out);
    }
    /*
    if (isFilename (met_dual_fn)) {
        // See if we have the right type of solve
        GrbLPSolver* lpsolver = dynamic_cast<GrbLPSolver*> (solver.get());
        if (lpsolver != nullptr) {
            cout << "Writing " << met_dual_fn << endl;
            ofstream out (met_dual_fn);
            auto biomass_react = scenario.biomass_react();
            auto duals = lpsolver->get_dual_vals();
            
            out << "# Produced by " << progname << " on " << todayString() << endl;
            out << "# datafile "  << data_fn << endl;
            out << "# " << summary.str() << endl;
            out << "# <met> <is-biomass-met> <lb-dual> <ub-dual>" << endl;
            for (auto& met_ptr : scenario.metabolites()) {
                auto met = met_ptr.get();
                bool is_biomass =
                    biomass_react->uses(met) &&
                    biomass_react->met_coeff (met) < 0;
                
                // Negate duals for human consumption
                out << met->name() << " " << is_biomass << " " <<
                    -duals->lb_dual(met) << " " << -duals->lb_dual(met) << endl;
            }
        }
    }
    */
    
    if (isFilename (cs_dot_fn)) {
        cout << "Writing " << cs_dot_fn << endl;
        ofstream dot (cs_dot_fn);
        dot << "/* Produced by " << progname << " on " << todayString() <<
            " */"<< endl;
        dot << "/* datafile "  << data_fn << " */" << endl;
        dot << "/* " << summary.str() << " */" << endl;
        
        if (sol != nullptr)
            sol->write_carbon_source_dot (dot);
    }
    if (isFilename (sd_dot_fn)) {
        cout << "Writing " << sd_dot_fn << endl;
        ofstream dot (sd_dot_fn);
        dot << "/* Produced by " << progname << " on " << todayString() <<
            " */"<< endl;
        dot << "/* datafile "  << data_fn << " */" << endl;
        dot << "/* " << summary.str() << " */" << endl;
        
        if (sol != nullptr)
            sol->write_supply_demand_dot (dot);
    }
    if (isFilename (summary_fn)) {
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
        sol->draw_reduced_cost (dig.get());
        dig->waitAndWipe();
        scenario.draw_reachability (dig.get(), sol.get());
    }

    return 0;
}
