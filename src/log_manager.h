#ifndef SPHLATCH_LOG_MANAGER_H
#define SPHLATCH_LOG_MANAGER_H

/*
 *  log_manager.h
 *
 *
 *  Created by Andreas Reufer on 23.06.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"
#include <unistd.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <boost/lexical_cast.hpp>

#ifdef SPHLATCH_PARALLEL
#include "communication_manager.h"
#endif

namespace sphlatch
{
class LogManager
{
public:
typedef LogManager self_type;
typedef LogManager& self_reference;
typedef LogManager* self_pointer;

static self_reference instance(void);
static void destroy(void);

#ifdef SPHLATCH_PARALLEL
typedef sphlatch::CommunicationManager CommManagerType;
CommManagerType& CommManager;
#endif

protected:
LogManager(void);
~LogManager(void);
std::fstream locLog;

private:
double getAbsTime();
double getRelTime();

double startTime, relStartTime;

public:
void operator<<(std::string _str);
void zeroRelTime();
std::ostringstream stream;
void flushStream();


static self_pointer _instance;
};

LogManager::self_pointer LogManager::_instance = NULL;

LogManager::self_reference LogManager::instance(void)
{
  if (_instance != NULL)
    return *_instance;
  else
    {
      _instance = new LogManager;
      return *_instance;
    }
}

///
/// here we need an explicit destructor, in order to properly close
/// the log files. calling the standard destructor doens't really work
/// in the singleton paradigm
///
void LogManager::destroy()
{
  delete _instance;
  _instance = NULL;
}

LogManager::LogManager(void)
#ifdef SPHLATCH_PARALLEL
  : CommManager(CommManagerType::instance())
#endif
{
#ifdef SPHLATCH_PARALLEL
  const size_t myDomain = CommManager.getMyDomain();
#else
  const size_t myDomain = 0;
#endif
  std::string logFilename = "logDomain000";
  std::string domString = boost::lexical_cast<std::string>(myDomain);
  logFilename.replace(logFilename.size() - 0 - domString.size(),
                      domString.size(), domString);
  locLog.open(logFilename.c_str(), std::ios::out);

  /// determine start time
  relStartTime = 0;
  startTime = 0;
  startTime = getAbsTime();

  /// get the hostname (unportable syscall,
  /// does not work on non-POSIX systems)
  char hostnameBuff[1024];
  gethostname(hostnameBuff, 1024);

  std::ostringstream header;
  header << "start log for domain " << myDomain << " on " << hostnameBuff << " (pid " << getpid() << ")";

  (*this) << header.str();
}

LogManager::~LogManager(void)
{
  (*this) << "close log";
  locLog.close();
}

void LogManager::operator<<(std::string _str)
{
  locLog << std::fixed << std::right << std::setw(18) << std::setprecision(6)
                         << getRelTime() << "      " << _str << "\n" << std::flush;
  relStartTime = getAbsTime();
}

void LogManager::zeroRelTime()
{
  relStartTime = getAbsTime();

  locLog << std::fixed << std::right << std::setw(18) << std::setprecision(6)
         << getAbsTime() << "***\n" << std::flush;
}


void LogManager::flushStream()
{
  stream << std::flush;
  (*this) << stream.str();
  stream.str("");
}

double LogManager::getRelTime()
{
  return(getAbsTime() - relStartTime);
}

double LogManager::getAbsTime()
{
#ifdef SPHLATCH_PARALLEL
  return CommManager.wtime() - startTime;
#else
  return 0.; /// todo: use another timer
             /// for serial version;
#endif
}
};

#endif
