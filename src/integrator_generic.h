#ifndef BHTREE_INTEGRATOR_GENERIC_H
#define BHTREE_INTEGRATOR_GENERIC_H

/*
 *  integrator_generic.h
 *
 *  A MetaIntegrator needs to have the following functions:
 *   bootstrap():   used to bootstrap the integrator, the derivatives should
 *                  contain the correct values, so that they can be saved
 *   integrate(dt): the integration function advancing the system with dt. the
 *                  derivatives should not be zeroed and contain the correct
 *                  values
 *
 *  Created by Andreas Reufer on 26.05.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

//#include "communication_manager.h"
#include "typedefs.h"
#include "particle_manager.h"
#ifdef SPHLATCH_PARALLEL
#include "communication_manager.h"
#endif

namespace sphlatch {
///
/// the integrator base class
///
/// virtual functions are usually considered slow, but
/// that's not so bad here
///
class GenericIntegrator {
typedef sphlatch::ParticleManager partManagerType;

#ifdef SPHLATCH_PARALLEL
typedef sphlatch::CommunicationManager commManagerType;
#endif

public:
GenericIntegrator(void) :
  PartManager(partManagerType::instance())
#ifdef SPHLATCH_PARALLEL
  ,CommManager(commManagerType::instance())
#endif
{
};
virtual ~GenericIntegrator()
{};

public:

protected:
partManagerType& PartManager;
#ifdef SPHLATCH_PARALLEL
commManagerType& CommManager;
#endif
};


class GenericMetaIntegrator {
typedef sphlatch::ParticleManager partManagerType;

public:
GenericMetaIntegrator() :
  PartManager(partManagerType::instance())
{};

virtual ~GenericMetaIntegrator(void)
{};

protected:
partManagerType& PartManager;

};
};

#endif
