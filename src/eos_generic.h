#ifndef SPHLATCH_EOS_GENERIC
#define SPHLATCH_EOS_GENERIC

/*
 *  eos_generic.h
 *
 * 
 *  Created by Andreas Reufer on 26.07.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"
#include "particle_manager.h"
#ifdef SPHLATCH_LOGGER
#include "log_manager.h"
#endif

namespace sphlatch
{

class EOS
{

protected:
typedef sphlatch::ParticleManager partManagerType;
partManagerType& PartManager;
#ifdef SPHLATCH_LOGGER
typedef sphlatch::LogManager logManagerType;
logManagerType& Logger;
#endif

valvectRefType rho, p, u;
idvectRefType mat;

public:
EOS() :
PartManager(partManagerType::instance()),
#ifdef SPHLATCH_LOGGER
Logger(logManagerType::instance()),
#endif
rho(PartManager.rho),
p(PartManager.p),
u(PartManager.u),
mat(PartManager.mat)
{};

~EOS() {};

};

}
#endif

