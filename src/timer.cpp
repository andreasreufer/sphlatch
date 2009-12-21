#ifndef TIMER_CPP
#define TIMER_CPP

/*
 *  timer.cpp
 *
 *
 *  Created by Andreas Reufer on 12.04.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#ifdef SPHLATCH_OPENMP
 #include <omp.h>
#elif SPHLATCH_MPI
 #include <mpi.h>
#else
 #include <sys/time.h>
#endif

namespace sphlatch {
class Timer {
public:
   Timer() { }
   ~Timer() { }

private:
   double startTime, roundTime;

public:
   double start()
   {
      const double curTime = getTime();

      startTime = curTime;
      roundTime = curTime;
      return(startTime);
   }

   double lapse()
   {
      const double curTime = getTime();

      roundTime = curTime;
      return(curTime - startTime);
   }

   double getRoundTime()
   {
      const double curTime = getTime();
      const double relTime = curTime - roundTime;

      roundTime = curTime;
      return(relTime);
   }

private:
   double getTime()
   {
#ifdef SPHLATCH_OPENMP
      return(omp_get_wtime());

#elif SPHLATCH_MPI
      return(MPI::Wtime());

#else
      struct timeval  tv;
      struct timezone tz;
      gettimeofday(&tv, &tz);
      return(static_cast<double>(tv.tv_sec) + static_cast<double>(tv.tv_usec) /
             1.e6);
#endif

   }
};
};

#endif
