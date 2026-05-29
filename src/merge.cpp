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
    string scen1_fn = "";
    string scen2_fn = "";
    string out_fn = "out.pld";
    string debug_str = "";

    // Keep scores from which scenario
    int which_scen = 1;

    Params params;

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

    Scenario scen1;
    Scenario scen2;
    
    cout << "Reading data from " << scen1_fn << endl;
    scen1.read_data (scen1_fn);

    cout << "Reading data from " << scen2_fn << endl;
    scen2.read_data (scen2_fn);

    cout << "Compare" << endl;

    // Translation from scen1 metabolite to scen2
    vector<const Metabolite*> xlate_met_12 (scen1.num_metabolites());
    // Translation from scen2 metabolite to scen1
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
