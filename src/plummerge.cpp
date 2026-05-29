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

// Return a "standardised" name from a string
string std_name (string str)
{
    string name = toLower (str);
    for (auto it = name.begin(); it != name.end(); it++) {
        switch (*it) {
        case '[': *it = '('; break;
        case ']': *it = ')'; break;
        case '_': *it = '-'; break;
        }
    }
    return name;
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
    
    // Ignore differences in metabolite full name
    bool ignore_met_descr = false;
    // exclude ex/dm from scen1?
    bool exclude_ex_dm_1 = false;
    // include ex/dm from scen2?
    bool exclude_ex_dm_2 = false;

    Params params;

    Opts opts (R"EOF(
    Merge two plum scenarios.
    - If metabolite full name is different, different versions of the met are
      created with unique-ified names. This can be switch off with '-md'
    - Reactions with the same name are checked for equality
    - If both scenarios contain a metabolite or reaction,
      '-w' switch controls which version is used
    - In particular '-w' switch controls reaction costs.
)EOF"
    );

    opts.add_opt ("-w", &which_scen, "Use scores/data from which scenario");
    opts.add_opt ("-md", &ignore_met_descr, "Ignore differences in metabolite description");
    opts.add_opt (
        "-ex1", &exclude_ex_dm_1, "Exclude EX/DM reactions from scenario 1?"
    );
    opts.add_opt (
        "-ex2", &exclude_ex_dm_2, "Exclude EX/DM reactions from scenario 2?"
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

    if (which_scen != 1 && which_scen != 2)
        opts.usage ("-w must be 1 or 2");

    Scenario scen1;
    Scenario scen2;
    
    cout << "Reading data from " << scen1_fn << endl;
    scen1.read_data (scen1_fn);

    cout << "Reading data from " << scen2_fn << endl;
    scen2.read_data (scen2_fn);

    cout << "Compare" << endl;

    string scen1_uniq = "_" + basename(scen1_fn).substr (0,2);
    string scen2_uniq = "_" + basename(scen2_fn).substr (0,2);
    if (equal_ic(scen1_uniq,scen2_uniq)) {
        scen1_uniq += "1";
        scen2_uniq += "2";
    }

    list<string> error;
    int error_count = 0;
    list<string> warn;
    int warn_count = 0;
    int prob_met_count = 0;
    int prob_react_count = 0;

    int excluded_ex_dm_count = 0;
    
    // Indicators of problematicb metabolites (differs in descr)
    vector<bool> prob_met1 (scen1.num_metabolites(), false);
    vector<bool> prob_met2 (scen2.num_metabolites(), false);

    // Translation from scen1 metabolite to scen2
    vector<const Metabolite*> xlate_met_12 (scen1.num_metabolites(), nullptr);
    // Translation from scen2 metabolite to scen1
    vector<const Metabolite*> xlate_met_21 (scen2.num_metabolites(), nullptr);

    // Match in 2 for reaction in 1
    vector<const Reaction*> match_react1 (scen1.num_reactions(), nullptr);
    // Match in 1 for reaction in 2
    vector<const Reaction*> match_react2 (scen2.num_reactions(), nullptr);

    for (auto& met1 : scen1.metabolites()) {
        DEBUG ('a', "Met1 " << met1->name());
        auto met2 = (Metabolite*) scen2.find_metabolite (met1->name());
        xlate_met_12[met1->index()] = met2;
        if (met2 == nullptr) {
            DEBUG ('b', "  No match");
            continue;
        }
        xlate_met_21[met2->index()] = met1.get();
        DEBUG ('b', "  Matched met2 " << met2->name());
        string name1 = std_name (met1->full_name());
        string name2 = std_name (met2->full_name());
        size_t len1 = name1.length();
        size_t len2 = name2.length();
        size_t minlen = std::min(len1,len2);
        if (
            name1.substr(0,minlen).compare(name2.substr(0,minlen)) != 0
        ) {
            DEBUG ('b', "    Descr differs");
            DEBUG ('b', "      minlen " << minlen);
            DEBUG ('b', "      met1: " << met1->full_name() << " " << name1);
            DEBUG ('b', "      met2: " << met2->full_name() << " " << name2);
            warn_count++;
            prob_met_count++;
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
            );
            if (ignore_met_descr) {
                if (which_scen == 2 || (len2 > len1)) {
                    // Replace name used in scen 1
                    met1->set_name (met2->name());
                }
            }
            else {
                prob_met1[met1->index()] = true;
                prob_met2[met2->index()] = true;
            }
        }
    }
    
    vector<bool> met_used1 (scen1.num_metabolites(), false);
    vector<bool> met_used2 (scen2.num_metabolites(), false);
    
    for (auto& react1 : scen1.reactions()) {
        DEBUG ('a', "React1 " << react1->name());

        if (exclude_ex_dm_1 && react1->is_ex_dm()) {
            DEBUG ('a', "  Is EX/DM - skip");
            continue;
        }
        
        auto react2 = (Reaction*)scen2.find_reaction (react1->name());

        if (react2 == nullptr) {
            DEBUG ('b', "  No match");
            continue; 
        }
        if (exclude_ex_dm_2 && react2->is_ex_dm()) {
            DEBUG ('a', "  React2 is EX/DM - skip");
            continue;
        }
        DEBUG ('b', "  Matched " << react2->name());
        
        match_react1[react1->index()] = react2;
        match_react2[react2->index()] = react1.get();
        
        bool ok = true;
        string error_met = "";
        // Does this reaction reference a problematic metabolite?
        bool prob = false;
        string prob_met = "";
        
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
                else if (prob_met1[idx1]) {
                    prob = true;
                    prob_met = met1->name();
                    match_react1[react1->index()] = nullptr;
                    match_react2[react2->index()] = nullptr;
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
            DEBUG ('a', "    Problem with " << error_met);
            DEBUG ('a', "      filenames " << scen1_fn << " " << scen2_fn);
            DEBUG ('a', "      react1 formula " << react1->formula());
            DEBUG ('a', "      react2 formula " << react2->formula());

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
            
            react1->set_name (react1->name() + scen1_uniq);
            react2->set_name (react2->name() + scen2_uniq);
        }
        else if (prob) {
            DEBUG ('a', "    Problematic met " << prob_met);
            DEBUG ('a', "      filenames " << scen1_fn << " " << scen2_fn);
            DEBUG ('a', "      react1 formula " << react1->formula());
            DEBUG ('a', "      react2 formula " << react2->formula());

            prob_react_count++;
            warn_count++;
            warn.push_back (
                "Reaction " + react1->name() +
                " references problematic metabolite " + prob_met
            );
            warn.push_back (
                limeFormat (
                    "%20s: %s", scen1_fn.c_str(), react1->formula().c_str()
                )
            );
            warn.push_back (
                limeFormat (
                    "%20s: %s", scen2_fn.c_str(), react2->formula().c_str()
                )
            );

            react1->set_name (react1->name() + scen1_uniq);
            react2->set_name (react2->name() + scen2_uniq);
        }
    }

    DEBUG ('a', "Check reactions in scen2");
    for (auto& react2 : scen2.reactions()) {
        DEBUG (
            'a', react2->name() << " num in " << react2->in_mets().size() <<
            " num out " << react2->out_mets().size()
        );
        if (exclude_ex_dm_2 && react2->is_ex_dm()) {
            DEBUG ('a', react2->name() << " is EX/DM - skip");
            continue;
        }
        
        for (auto idx2 : react2->mets()) {
            met_used2[idx2] = true;
        }
    }

    // Write out merged scenario
    DEBUG ('a', "Write merged scenario");
    size_t num_reacts = 0;
    size_t num_mets = 0;
    ofstream out (out_fn);
    out << "# Produced from " << scen1_fn << " and " <<  scen2_fn <<
        " on " << todayString() << " by " << progname << endl;
    out << "# Favour scenario " << which_scen <<
        ", ignore metabolite description differences: " <<
        (ignore_met_descr ? "true" : "false") << endl;
    
    for (auto& met1 : scen1.metabolites()) {
        if (!met_used1[met1->index()])
            continue;

        if (prob_met1[met1->index()]) {
            met1->set_name (met1->name() + scen1_uniq);
            met1->write_to (out);
        }
        else if (
            xlate_met_12[met1->index()] == nullptr ||
            which_scen == 1
        ) {
            num_mets++;
            met1->write_to (out);
        }
        else {
            num_mets++;
            xlate_met_12[met1->index()]->write_to (out);
        }
    }
    // Write out mets in scen2 that have no equiv in scen1
    for (auto& met2 : scen2.metabolites()) {
        if (!met_used2[met2->index()])
            continue;
        if (prob_met2[met2->index()]) {
            num_mets++;
            met2->set_name (met2->name() + scen2_uniq);
            met2->write_to (out);
        }
        else {
            auto met1 = xlate_met_21[met2->index()];
            if (met1 == nullptr) {
                num_mets++;
                met2->write_to (out);
            }
            else {
                // Double check we were not relying on it
                if (!met_used1[met1->index()]) {
                    // Huh. Scen 1 uses it but not scen 2. Write it out
                    num_mets++;
                    met2->write_to (out);
                }
            }
        }
    }

    for (auto& react1 : scen1.reactions()) {
        if (exclude_ex_dm_1 && react1->is_ex_dm()) {
            excluded_ex_dm_count++;
            const char* exdm = "EX";
            if (react1->is_dm())
                exdm = "DM";
            warn.push_back (
                limeFormat (
                    "Reaction %30s from %20s is %s - skipped",
                    react1->name().c_str(), scen1_fn.c_str(), exdm
                )
            );
        }
        
        auto react2 = match_react1[react1->index()];
        
        if (
            react2 == nullptr ||
            which_scen == 1
        ) {
            num_reacts++;
            react1->write_to (out);
        }
    }
    // Write out reacts from scen2 that don't appear in scen1
    for (auto& react2 : scen2.reactions()) {
        if (exclude_ex_dm_2 && react2->is_ex_dm()) {
            excluded_ex_dm_count++;
            const char* exdm = "EX";
            if (react2->is_dm())
                exdm = "DM";
            warn.push_back (
                limeFormat (
                    "Reaction %30s from %20s is %s - skipped",
                    react2->name().c_str(), scen2_fn.c_str(), exdm
                )
            );
        }
        
        auto react1 = match_react2[react2->index()];
        
        if (
            react1 == nullptr ||
            which_scen == 2
        ) {
            num_reacts++;
            react2->write_to (out);
        }
    }

    string err_fn = "plummerge.out";
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
        cout << "including " << prob_met_count <<
            " problematic metabolites, and" << endl;
        cout << "          " << prob_react_count <<
            " problematic reactions" << endl;
        cout << "          " << excluded_ex_dm_count <<
            " excluded EX/DM reactions" << endl;
        
        errout << warn_count << " warnings" << endl;
        errout << "including " << prob_met_count <<
            " problematic metabolites, and" << endl;
        errout << "          " << prob_react_count <<
            " problematic reactions" << endl;
        errout << "          " << excluded_ex_dm_count <<
            " excluded EX/DM reactions" << endl;

        if (prob_met_count > 0) {
            string have = "have";
            if (ignore_met_descr) 
                have = "have NOT";
            cout << "Note: problematic metabolites " << have <<
                " been replaced" << endl << endl;
            errout << "Note: problematic metabolites " << have <<
                " been replaced" << endl << endl;
        }
        
        for (auto text : warn) {
            cout <<  text << endl;
            errout <<  text << endl;
        }
    }
    cout << endl;
    cout << "Wrote merged scenario to " << out_fn << endl;
    cout << "                    with " << num_mets << " metabolites" << endl;
    cout << "                     and " << num_reacts << " reactions" << endl;
    cout << "         Wrote errors to " << err_fn << endl;
    
    return 0;
}
