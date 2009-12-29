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

/*#ifdef SPHLATCH_PARALLEL
 #include "communication_manager.h"
 #endif*/

#include "timer.cpp"

namespace sphlatch
{
class Logger
{
public:
   typedef Logger            self_type;
   typedef Logger&           self_reference;
   typedef Logger*           self_pointer;

   static self_reference instance(void);
   static void           destroy(void);
   typedef sphlatch::Timer   timerT;

/*#ifdef SPHLATCH_PARALLEL
   typedef sphlatch::CommunicationManager CommManagerType;
   CommManagerType& CommManager;
 #endif*/

public:
   Logger(void);
   ~Logger(void);
   std::fstream locLog;

private:
   timerT Timer;

public:
   void operator<<(std::string _str);

   std::ostringstream stream;
   void finishStep();
   void finishStep(std::string _str);
   void flushStream();

   static self_pointer _instance;
};

Logger::self_pointer Logger::_instance = NULL;

Logger::self_reference Logger::instance(void)
{
   if (_instance != NULL)
      return(*_instance);
   else
   {
      _instance = new Logger;
      return(*_instance);
   }
}

///
/// here we need an explicit destructor, in order to properly close
/// the log files. calling the standard destructor doens't really work
/// in the singleton paradigm
///
void Logger::destroy()
{
   delete _instance;
   _instance = NULL;
}

Logger::Logger(void)

/*#ifdef SPHLATCH_PARALLEL
   : CommManager(CommManagerType::instance())
 #endif*/
{
   //FIXME: in parallel version
   const size_t myDomain    = 0;
   std::string  logFilename = "logDomain000";
   std::string  domString   = boost::lexical_cast<std::string>(myDomain);

   logFilename.replace(logFilename.size() - 0 - domString.size(),
                       domString.size(), domString);
   locLog.open(logFilename.c_str(), std::ios::out);

   Timer.start();

   /// get the hostname (unportable syscall,
   /// does not work on non-POSIX systems)
   char hostnameBuff[1024];
   gethostname(hostnameBuff, 1024);

   std::ostringstream header;
   header << "start log for domain " << myDomain << " on " << hostnameBuff <<
   " (pid " << getpid() << ")";

   (*this) << header.str();
}

Logger::~Logger(void)
{
   (*this) << "close log";
   locLog.close();
}

void Logger::operator<<(std::string _str)
{
   locLog << std::fixed << std::right << std::setw(18) << std::setprecision(6)
                     << Timer.getRoundTime() << "      " << _str << "\n" <<
   std::flush;
}

void Logger::finishStep(std::string _str)
{
   locLog << std::fixed << std::right << std::setw(18) << std::setprecision(6)
          << Timer.lapse() << "***   " << _str << "\n" << std::flush;
}

void Logger::finishStep()
{
   finishStep("");
}

void Logger::flushStream()
{
   stream << std::flush;
   (*this) << stream.str();
   stream.str("");
}
};

#endif
