
#include "lime/error.h"

#include "mosh/lpsutil.h"

using namespace std;
using namespace mosh;

void mosh::check_lps_return (unsigned char retval, std::string where)
{
    if (retval != 1) 
        limeCrash ("Return val " + to_string(retval) + " from " + where);
}
