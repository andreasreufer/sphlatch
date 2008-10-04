#ifndef TIMER_CC
#define TIMER_CC

#include <ctime>

#include "Def.cc"

class Timer {
private:
  bool    running;
  clock_t beg, end, sum;

public:
  Timer() { running = false; reset(); }
  
  void reset() {
    if (running) { std::cerr << "Reset running timer!\n"; }
    sum = 0;
    running = false;
  }
 
  void start() {
    beg = clock();
    if (running) { std::cerr << "Timer already running!\n"; quit(1); }
    running = true;
  }
 
  void stop()  {
    end = clock();
    if (!running) { std::cerr << "Timer not started!\n"; quit(1); }
    sum += (end - beg);
    running = false;
  }
 
  float result() const { return (float)sum/(float)CLOCKS_PER_SEC; }
};

#endif
