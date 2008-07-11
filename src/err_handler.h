#ifndef SPHLATCH_ERR_HANDLER_H
#define SPHLATCH_ERR_HANDLER_H

/*
 *  err_handler.h
 *
 *
 *  Created by Andreas Reufer on 10.07.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"
#include "log_manager.h"
#include "particle_manager.h"

namespace sphlatch
{
class GenericError
{
public:
typedef sphlatch::LogManager logType;
typedef sphlatch::ParticleManager partManagerType;

GenericError();
~GenericError();

protected:
logType& Logger;
partManagerType& PartManager;
};

GenericError::GenericError(void)
  : Logger(logType::instance()), PartManager(partManagerType::instance())
{
  Logger << "SPHLATCH!!! some error has occured ...";
};

GenericError::~GenericError(void)
{
};
};

#endif
