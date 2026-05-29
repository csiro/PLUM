#pragma once

#include <vector>
#include <memory>

#include "lime/dig.h"
#include "lime/timekeeper.h"
#include "lime/constants.h"
#include "lime/rand.h"
#include "lime/debug.h"
#include "lime/strutil.h"

#include "mosh/scenario.h"
#include "mosh/params.h"
#include "mosh/solution.h"
#include "mosh/multisol.h"
#include "mosh/reaction.h"
#include "mosh/metabolite.h"

namespace mosh
{
    class GapSolver 
    {
    public:
        GapSolver(
            Scenario* scenario, Params* params, int seed
        ) :
            scenario_(scenario),
            params_(params),
            quiet_(false),
            rand_(seed),
            status_(UNSOLVED),
            best_sol_fn_(""),
            timer_()
        {
            if (params_->time_limit > 0)
                timer_.setTimeLimit (params_->time_limit);
        }

        size_t num_metabolites() const {return scenario_->num_metabolites();}
        const Metabolite* metabolite(size_t k)
        {
            return scenario_->metabolite(k);
        }
        
        size_t num_reactions() const {return scenario_->num_reactions();}
        Reaction* reaction(size_t k) const {return scenario_->reaction(k);}
        
        size_t num_experiments() const {return scenario_->num_experiments();}
        Experiment* experiment(size_t k) const
        {
            return scenario_->experiment(k);
        }
        
        std::string best_sol_fn() const {return best_sol_fn_;}
        void set_best_sol_fn (std::string best_sol_fn) {
            best_sol_fn_ = best_sol_fn;
        }

        int grb_seed() {return rand_.uniform0n_1 (1000000);}
        int rand_seed() {return rand_.generateSeed();}

        virtual SolutionPtr solve() = 0;

        virtual std::string summary() {return "";}

        double elapsed_time_secs () const {return timer_.elapsedTimeSecs();}
        bool has_time_left() const {return timer_.hasTimeLeft();}
        double time_left_seconds() const {return timer_.timeLeftSecs();}
        StatusEnum status() const {return status_;}
        std::string status_str() const
        {
            std::vector<std::string> stat_str =
                {
                    "unsolved",
                    "cts_optimal",
                    "int_optimal",
                    "suboptimal",
                    "infeasible",
                    "random"
                };
            return stat_str[status_];
        }

        bool quiet() const {return quiet_;}
        virtual void set_quiet (bool quiet) {
            quiet_ = quiet;
        }
        void on_screen (std::string message) {
            DEBUG ('A', message);
            if (!quiet_) {
                std::cout <<
                    escSeqColour (lime::CYAN) <<
                    message <<
                    escSeqColour (lime::RESET) <<
                    std::endl;
            }
        }

    protected:
        Scenario* scenario_;
        Params* params_;
        bool quiet_;

        lime::Rand rand_;

        StatusEnum status_;
        
        std::string best_sol_fn_;
        
    private:
        lime::TimeKeeper timer_;
    };

    using GapSolverPtr = std::shared_ptr<GapSolver>;
}
