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
#include "io_manager.h"

namespace sphlatch
{
class GenericError
{
public:
typedef sphlatch::LogManager logType;
typedef sphlatch::ParticleManager partManagerType;
typedef sphlatch::IOManager ioManagerType;

GenericError();
~GenericError();

protected:
logType&         Logger;
partManagerType& PartManager;
ioManagerType&   IOManager;
};

GenericError::GenericError(void)
  : Logger(logType::instance()),
    PartManager(partManagerType::instance()),
    IOManager(ioManagerType::instance())
{
  Logger << "SPHLATCH!!! an error has occured ...";
};

GenericError::~GenericError(void)
{
};
};

#endif
