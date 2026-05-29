
#include "lime/error.h"

#include "mosh/highsutil.h"


using namespace std;
using namespace mosh;

void mosh::check_highs_status (HighsStatus status, std::string where)
{
    switch (status) {
    case HighsStatus::kOk:
        // Do nothing
        break;
    case HighsStatus::kWarning:
        limeWarning ("HiGHS warning from " + where);
        break;
    case HighsStatus::kError:
        limeWarning ("HiGHS error from " + where);
        break;
    }
}
