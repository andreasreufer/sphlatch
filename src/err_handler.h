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

#include <string>

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


class FileNotFound : public GenericError
{
public:
FileNotFound(std::string _filename);
~FileNotFound();

};

FileNotFound::FileNotFound(std::string _filename)
{
  Logger.stream << "file not found: " << _filename;
  Logger.flushStream();
};

FileNotFound::~FileNotFound()
{
};

};
#endif
