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
#include "particle_manager.h"
#include "io_manager.h"
#ifdef SPHLATCH_LOGGER
#include "log_manager.h"
#endif

namespace sphlatch
{
class GenericError
{
public:
#ifdef SPHLATCH_LOGGER
typedef sphlatch::LogManager logType;
#endif
typedef sphlatch::ParticleManager partManagerType;
typedef sphlatch::IOManager ioManagerType;

GenericError();
~GenericError();

protected:
#ifdef SPHLATCH_LOGGER
logType&         Logger;
#endif
partManagerType& PartManager;
ioManagerType&   IOManager;
};

GenericError::GenericError(void) :
#ifdef SPHLATCH_LOGGER
    Logger(logType::instance()),
#endif
    PartManager(partManagerType::instance()),
    IOManager(ioManagerType::instance())
{
#ifdef SPHLATCH_LOGGER
  Logger << "SPHLATCH!!! an error has occured ...";
#else
  std::cerr << "SPHLATCH!!! an error has occured ...\n";
#endif
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
#ifdef SPHLATCH_LOGGER
  Logger.stream << "file not found: " << _filename;
  Logger.flushStream();
#else
  std::cerr << "file not found: " << _filename << "\n";
#endif
};

FileNotFound::~FileNotFound()
{
};

};
#endif
