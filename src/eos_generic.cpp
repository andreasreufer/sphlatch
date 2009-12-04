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
#ifdef SPHLATCH_LOGGER
#include "logger.cpp"
#endif

namespace sphlatch
{

class EOS
{
protected:
#ifdef SPHLATCH_LOGGER
typedef sphlatch::Logger logT;
logT& Logger;
#endif

public:
EOS()
#ifdef SPHLATCH_LOGGER
: Logger(logT::instance()),
#endif
{};

~EOS() {};

};

}
#endif

