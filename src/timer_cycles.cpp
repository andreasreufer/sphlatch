#ifndef TIMER_CPP
#define TIMER_CPP

/*
 *  timer_cycles.cpp
 *
 *
 *  Created by Andreas Reufer on 12.04.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

namespace sphlatch {
class CycleTimer {
public:
   CycleTimer() { }

   ~CycleTimer() { }

private:
   long long unsigned int        startCycles;
   static long long unsigned int cyclBuff;

public:
   void start()
   {
      getCounts();
      startCycles = cyclBuff;
   }

   long long unsigned int lapse()
   {
      getCounts();
      return(cyclBuff - startCycles);
   }

private:
   void getCounts()
   {
      // read time stamp counter on x86 or amd64 CPU
      // warning: not portable!
      __asm volatile ("rdtsc" : "=A" (cyclBuff));
   }
};

long long unsigned int CycleTimer::cyclBuff;
};

#endif
