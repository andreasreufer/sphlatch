#ifndef SPHLATCH_EOS_GENERIC
#define SPHLATCH_EOS_GENERIC

/*
 *  eos_generic.h
 *
 * 
 *  Created by Andreas Reufer on  3.12.09 
 *  Copyright 2009 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"
//#ifdef SPHLATCH_LOGGER
//#include "log_manager.h"
//#endif

namespace sphlatch
{

class EOS
{

protected:
/*#ifdef SPHLATCH_LOGGER
typedef sphlatch::LogManager logManagerType;
logManagerType& Logger;
#endif*/


public:
EOS() :
/*#ifdef SPHLATCH_LOGGER
Logger(logManagerType::instance()),
#endif*/
{};

~EOS() {};

};

}
#endif

