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
#include "log_manager.h"


namespace sphlatch
{

class EOS
{

protected:
typedef sphlatch::ParticleManager partManagerType;
partManagerType& PartManager;
typedef sphlatch::LogManager logManagerType;
logManagerType& Logger;

valvectRefType rho, p, u;
idvectRefType mat;

public:
EOS() :
PartManager(partManagerType::instance()),
Logger(logManagerType::instance()),
rho(PartManager.rho),
p(PartManager.p),
u(PartManager.u),
mat(PartManager.mat)
{};

~EOS() {};

};

}
#endif

