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
#include "lime/sortpairt.h"

#include "mosh/builddate.h" // For _build_date
#include "mosh/gapsolver.h" 
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

SolutionPtr
read_file (Scenario* scenario, string fn, vector<int>& count)
{
    Params params;
    SolutionPtr sol = make_shared<Solution> (scenario, &params);
    
    LineReader reader (fn);

    string line;
    while (reader.getLine (line)) {
        LimeTok tok (line);
        bool error = false;
        string name = tok.nextString(error);
        double flux = tok.nextDouble(error);

        auto react = scenario->find_reaction (name);
        if (react == nullptr)
            reader.error (string("Lost reaction ") + name);
        sol->set_flux (react, flux);
        count[react->index()] += 1;
    }
    return sol;
}


int
main (int argc, const char* argv[]) 
{
    string data_fn = "";
    string dig_fn = "overlap.dig";
    vector<string> sol_fn;
    string debug_str = "";

    Params params;

    Opts opts (R"EOF(
    Plot solutions to gap-filling problem
)EOF"
    );

    opts.add_opt_filename ("-g", &dig_fn, "Dig filename");
    opts.add_opt ("-d", &debug_str, "Debug string");
    opts.add_opt ("-j", &params.init_biomass_obj_mult, "Biomass objective multiplier");

    opts.add_arg ("data.pld", &data_fn, "Input data");
    
    opts.add_arg ("sol.out...", &sol_fn, "Solution file(s)");
    

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

    Scenario scenario;
    
    cout << "Reading data from " << data_fn << endl;
    scenario.read_data (data_fn);
    //scenario.add_dummy_reactions();

    vector<int> count (scenario.num_reactions(), 0);
    vector<SolutionPtr> sols;

    for (auto fn : sol_fn) {
        cout << "Read sol from " << fn << endl;

        auto sol = read_file (&scenario, fn, count);
        sols.push_back (sol);
    }
    double max_flux = 0.0f;           
    vector<double> max_react_flux (scenario.num_reactions(), 0.0f);
    for (auto& sol : sols)  {
        for (size_t k = 0; k < scenario.num_reactions(); k++) {
            if (sol->flux(k) > max_flux)
                max_flux = sol->flux(k);
            if (sol->flux(k) > max_react_flux[k])
                max_react_flux[k] = sol->flux(k);
        }
    }

    vector<size_t> proportion (sols.size() + 1, 0);
    size_t num_used = 0;
    for (size_t k = 0; k < scenario.num_reactions(); k++) {
        if (!scenario.reaction(k)->is_dummy()) {
            proportion[count[k]] += 1;
            if (count[k] > 0)
                num_used++;
        }
    }
    cout << "Proportions" << endl;
    cout << "Count   of used   overall " << endl;
//          "123 2345678901 234567890"
    for (size_t i = 0; i <= sols.size(); i++) {
        int ofused =
            (int) round (100.0f * proportion[i] / num_used);
        int overall =
            (int) round (
                100.0f * proportion[i] / scenario.orig_num_reactions()
            );
        if (i == 0)
            ofused = 0;
        cout << setw(3) << i << 
            setw(11) << ofused << '%' << 
            setw(9) << overall << '%' << endl;
    }

    cout << "Writing to " << dig_fn << endl;
    Dig dig (dig_fn);
    dig.showMessage (
        string("Produced by ") + progname + " on " + todayString()
    );
    dig.showMessage (string ("datafile ") + data_fn);

    dig.title ("Solution overlap");
    dig.setPrecision (8);

    for (int iter = 0; iter < 3; iter++) {
        // Iter 0: y axis is simply sol number
        // Iter 1: y axis is relative obj
        // Iter 2: y axis is absolute obj
        switch (iter) {
        case 1: dig.yLabel ("Relative Objective"); break;
        case 2: dig.yLabel ("Absolute Objective"); break;
        }

        vector<double> xtic (scenario.num_reactions(), 0.0f);
        int num = 0;
        for (size_t k = 0; k < scenario.num_reactions(); k++) {
            if (count[k] > 0) {
                xtic[k] = (double)num;
                dig.xTic (xtic[k], scenario.reaction(k)->name());
                num++;
            }
        }
        if (iter == 0) {
            // Y tics
            //for (size_t i = 0; i < sol_fn.size(); i++) 
                //dig.yTic (i+1, sol_fn[i]);
        }
    
        double base_diam = 30.0f;
        for (size_t i = 0; i < sols.size(); i++) {
            auto& sol = sols[i];
            dig.style(i+1);
            dig.label (sol_fn[i]);
            double rel_obj = sol->rel_obj_value();
            double abs_obj = sol->abs_obj_value();
            
            double y = 0;
            switch (iter) {
            case 0: y = (double)(i+1); break;
            case 1: y = sol->rel_obj_value(); break;
            case 2: y = sol->abs_obj_value(); break;
            }
            for (size_t k = 0; k < scenario.num_reactions(); k++) {
                if (sol->flux(k) > 0.0f) {
                    double x = xtic[k];
                    int diam = base_diam * sqrt(sol->flux(k)) / sqrt(max_flux);
                    dig.circle (x, y, diam, true);
                }
            }
        }
        dig.waitAndWipe();
    }
    
    dig.yLabel ("Count"); 
    int num = 0;
    vector<double> xtic (scenario.num_reactions(), 0.0f);
    for (size_t k = 0; k < scenario.num_reactions(); k++) {
        if (count[k] > 0) {
            xtic[k] = (double)num;
            dig.xTic (xtic[k], scenario.reaction(k)->name());
            num++;
        }
    }
    vector<int> count2 (scenario.num_reactions(), 0);
    for (size_t i = 0; i < sols.size(); i++) {
        dig.style(i+1);
        dig.label (sol_fn[i]);
        for (size_t k = 0; k < scenario.num_reactions(); k++) {
            if (sols[i]->flux(k) > 0) {
                double x = xtic[k];
                double y1 = (double)count2[k];
                double y2 = y1 + 1;
                dig.draw (x, y1, x, y2);
                count2[k]++;
            }
        }
    }

    dig.waitAndWipe();
    dig.yLabel ("Flux"); 
    for (size_t k = 0; k < scenario.num_reactions(); k++) {
        if (count[k] > 0) {
            dig.xTic (xtic[k], scenario.reaction(k)->name());
        }
    }
    for (size_t i = 0; i < sols.size(); i++) {
        dig.style(i + 1, Dig::Markers::HBAR);
        dig.style(i+1);
        dig.label (sol_fn[i]);
    }
    double step = 1.0f / (sols.size() + 1);
    for (size_t k = 0; k < scenario.num_reactions(); k++) {
        if (count[k] > 0) {
            double x = xtic[k];
            for (size_t i = 0; i < sols.size(); i++) {
                double y = sols[i]->flux(k);
                dig.style(i + 1, Dig::Markers::HBAR);
                dig.draw (x, 0, x, y);
                x += step;
            }
        }
    }
    dig.wait();
    SortPairT<size_t> sorter;
    for (size_t k = 0; k < scenario.num_reactions(); k++) {
        if (count[k] > 0) {
            double x = xtic[k];
            double max_diff = 0.0f;
            for (size_t i = 0; i < sols.size(); i++) {
                double y = sols[i]->flux(k);
                if (sols[i]->flux(k) < max_react_flux[k]) {
                    dig.style(
                        sols.size() + 1,
                        Dig::Markers::HBAR,
                        Dig::StyleStroke::THICK
                    );
                    double y2 = max_react_flux[k];
                    dig.draw (x, y, x, y2);

                    double diff = max_react_flux[k] - sols[i]->flux(k);
                    if (diff > max_diff)
                        max_diff = diff;
                }                    
                x += step;
            }
            if (max_diff > 0.0f)
                sorter.add (k, -max_diff);
        }
    }

    // Redo - sorted by difference
    dig.waitAndWipe();
    dig.yLabel ("Flux");
    for (size_t j = 0; j < sorter.size(); j++) {
        auto k = sorter[j];
        dig.xTic (j, scenario.reaction(k)->name());
    }
    for (size_t i = 0; i < sols.size(); i++) {
        dig.style(i + 1, Dig::Markers::HBAR);
        dig.style(i+1);
        dig.label (sol_fn[i]);
    }
    for (size_t j = 0; j < sorter.size(); j++) {
        auto k = sorter[j];
        double x = j;
        for (size_t i = 0; i < sols.size(); i++) {
            double y = sols[i]->flux(k);
            dig.style(i + 1, Dig::Markers::HBAR);
            dig.draw (x, 0, x, y);
            x += step;
        }
    }
    dig.wait();
    for (size_t j = 0; j < sorter.size(); j++) {
        auto k = sorter[j];
        double x = j;
        for (size_t i = 0; i < sols.size(); i++) {
            double y = sols[i]->flux(k);
            if (sols[i]->flux(k) < max_react_flux[k]) {
                dig.style(
                    sols.size() + 1,
                    Dig::Markers::HBAR,
                    Dig::StyleStroke::THICK
                );
                double y2 = max_react_flux[k];
                dig.draw (x, y, x, y2);
            }                    
            x += step;
        }
    }

    // Calculate the max difference between first and rest of solutions,
    // for use as ordering
    vector<double> diff (scenario.num_metabolites());
    vector<double> metbal0 (scenario.num_metabolites());
    sols[0]->calc_metbal (metbal0);
    vector<double> metbal (scenario.num_metabolites());
    for (size_t i = 1; i < sols.size(); i++) {
        sols[i]->calc_metbal (metbal);
        for (size_t k = 0; k < scenario.num_metabolites(); k++) {
            double d = fabs(metbal[k] - metbal0[k]);
            if (d > diff[k])
                diff[k] = d;
        }
    }
    sorter.clear();
    for (size_t k = 0; k < scenario.num_metabolites(); k++) {
        // Negate diff to make descending sort
        sorter.add (k, -diff[k]);
    }

    dig.waitAndWipe();
    dig.title ("Metabolite Balance"); 
    dig.yLabel ("Supply/Residual"); 
    for (size_t k = 0; k < scenario.num_metabolites(); k++) {
        size_t met_k = sorter[k];
        dig.xTic ((double)(k+1), scenario.metabolite(met_k)->name());
    }
    for (size_t i = 0; i < sols.size(); i++) {
        dig.style(i + 1, Dig::Markers::HBAR);
        dig.style(i+1);
        dig.label (sol_fn[i]);
    }
    double offset = 0.0f;
    for (size_t i = 0; i < sols.size(); i++) {
        sols[i]->calc_metbal (metbal);
        for (size_t k = 0; k < scenario.num_metabolites(); k++) {
            size_t met_k = sorter[k];
            double x = (double)(k+1) + offset;
            double y = -metbal[met_k];
            dig.style(i + 1, Dig::Markers::HBAR);
            dig.draw (x, 0, x, y);
        }
        offset += step;
    }

    // Calculate the max difference between first and rest of solutions,
    // for use as ordering
    fill (diff.begin(), diff.end(), 0.0f);
    vector<double> metprod0 (scenario.num_metabolites());
    sols[0]->calc_metprod (metprod0);
    vector<double> metprod (scenario.num_metabolites());
    for (size_t i = 1; i < sols.size(); i++) {
        sols[i]->calc_metprod (metprod);
        for (size_t k = 0; k < scenario.num_metabolites(); k++) {
            double d = fabs(metprod[k] - metprod0[k]);
            if (d > diff[k])
                diff[k] = d;
        }
    }
    sorter.clear();
    for (size_t k = 0; k < scenario.num_metabolites(); k++) {
        // Negate diff to make descending sort
        sorter.add (k,-diff[k]);
    }

    dig.waitAndWipe();
    dig.title ("Metabolite Flux"); 
    dig.yLabel ("Production"); 
    for (size_t k = 0; k < scenario.num_metabolites(); k++) {
        size_t met_k = sorter[k];
        dig.xTic ((double)(k+1), scenario.metabolite(met_k)->name());
    }
    for (size_t i = 0; i < sols.size(); i++) {
        dig.style(i + 1, Dig::Markers::HBAR);
        dig.style(i+1);
        dig.label (sol_fn[i]);
    }
    offset = 0.0f;
    for (size_t i = 0; i < sols.size(); i++) {
        sols[i]->calc_metprod (metprod);
        for (size_t k = 0; k < scenario.num_metabolites(); k++) {
            size_t met_k = sorter[k];
            double x = (double)(k+1) + offset;
            double y = metprod[met_k];
            dig.style(i + 1, Dig::Markers::HBAR);
            dig.draw (x, 0, x, y);
        }
        offset += step;
    }
    
    
    return 0;
}
