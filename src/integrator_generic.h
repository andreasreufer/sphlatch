#ifndef BHTREE_INTEGRATOR_GENERIC_H
#define BHTREE_INTEGRATOR_GENERIC_H

/*
 *  integrator_generic.h
 *
 *
 *  Created by Andreas Reufer on 26.05.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

//#include "communication_manager.h"
#include "typedefs.h"
#include "particle_manager.h"
#include "communication_manager.h"

namespace sphlatch {
///
/// the integrator base class
///
/// virtual functions are usually considered slow, but
/// that's not so bad here
///
class GenericIntegrator {
typedef sphlatch::ParticleManager partManagerType;
typedef sphlatch::CommunicationManager commManagerType;

public:
GenericIntegrator(void) :
  PartManager(partManagerType::instance()),
  CommManager(commManagerType::instance())
{
};
virtual ~GenericIntegrator()
{};

public:

protected:
partManagerType& PartManager;
commManagerType& CommManager;
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
