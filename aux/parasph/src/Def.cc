#ifndef DEF_CC
#define DEF_CC

#include <iostream>
#include <cmath>
#include <string>
#include <unistd.h>
#include <stdlib.h>

#ifdef MPI_OK
#include <mpi.h>
#endif

#include "CommStrat.cc"

// define global variables
namespace global {
  bool      noisy;
  float     theta;
  int       dim, dim2, dim3, npro, rank, slice, totNumPart;
  CommStrat comm;

  void init() { comm.init(rank, npro); }

  void setTotNumPart(const int &n) {
    totNumPart   = n;
    slice        = (int)pow((pow(totNumPart, 1./3.) + 1), 2);
    dim          = (totNumPart - 1) / slice + 1;
    dim2         = dim * dim;
    dim3         = dim * dim2;
  }
}

// SPH data types
#ifdef DOUBLEPREC
#define ftype     double
#define MPI_ftype MPI::DOUBLE
#else
#define ftype     float
#define MPI_ftype MPI::FLOAT
#endif

// MPI shortcut (just for lazy people like me)
#define MCW MPI::COMM_WORLD

// Why isn't there a maximum function in standard C?
template <class T> 
T Max(const T &a, const T &b) { return a > b ? a : b; }
template <class T> 
T Max(const T &a, const T &b,
      const T &c)             { return Max(Max(a, b), c); }
template <class T> 
T Max(const T &a, const T &b,
      const T &c, const T &d) { return Max(Max(a, b), Max(c, d)); }
//template <class T>
//T Max(T *a, const int &size) {
//  T h = a[0];
//  for (int i = 1; i < size; i++) if (a[i] > h) h = a[i];
//  return h;
//}
template <class T> 
T Max(const T &a, const T &min, const std::string &warning) { 
  if (a < min) { std::cout << warning << " " << a << std::endl; return min; }
  return a;
}

// Instead of the maximum this function gives you its position
template <class T>
int maxLoc(T *a, const int &size) {
  int i, pos;
  for (i = 1, pos = 0; i < size; i++) if (a[i] > a[pos]) pos = i;
  return pos;
}

template <class T>
T Min(const T &a, const T &b) { return a < b ? a : b; }

template <class T>
T Fabs(const T &a) { return (T)fabs(a); }

// ... extending the String class ...
std::string trim(const std::string &input) {
  std::string str = input;
  str.erase(0, str.find_first_not_of(" "));
  str.erase(str.find_last_not_of(" ")+1); 
  return str;
}
std::string toUpper(const std::string &input) {
  std::string str = input;
  for (int i = 0; i < str.length(); i++) str[i] = toupper(str[i]);
  return str;
}

// Quit the simulation always with this function
void quit(int code) {
  char s [] = "Simulation stopped abnormally!\n";
  if (code) write(2, s, sizeof(s));
#ifdef MPI_OK
  if (MPI::Is_initialized()) MPI::Finalize();
#endif
  exit(code);
}

void listSwitches(const bool &noisy) {
#ifdef SWITCHES
  if (global::noisy) 
    std::cout << "compiler switches: " << SWITCHES << std::endl;
#else
  if (global::noisy) std::cout << "Warning: compiler variable SWITCHES not "
			       << "defined. It should contain a list of used "
			       << "switches." << std::endl;
#endif
}

// Check if all processors execute the same program :-)
#ifdef MPI_OK
#include "numberfile"
class ErrorVaryingBuildNumbers {};
void checkbuild(const int &rank, const bool &noisy) {
  int bn, res;
  if (rank == 0) bn = BUILDNUM;
  MCW.Bcast(&bn, 1, MPI::INT, 0);
  if (bn != BUILDNUM) bn = -1;
  MCW.Allreduce(&bn, &res, 1, MPI::INT, MPI_MIN);
  if (res < 0) throw ErrorVaryingBuildNumbers();
  else if (noisy) std::cout << "build #" << bn << std::endl;
}
#endif

// This is the timestep-condition-datatype which is used to administrate
// timesteps.
class dtcond {
public:
  double dt;
  ftype  val, dval;
  long   cond, num;
  enum { numcond = 12 };
  static const std::string condname[numcond];

  dtcond() { set(1.e30, -1, -1, -1.,-1.); }

  void set(const ftype &_dt,  const long  &_cond, 
	   const long  &_num, const ftype &_val, const ftype &_dval) {
    dt = _dt; cond = _cond; num = _num; val = _val; dval = _dval;
  }
  
  void print(const int &rank) {
    std::cout << "  timestep condition (" << condname[cond] << ") ";
    std::cout << rank << " " << dt << " ";
    if (cond > 0) std::cout << num << " value " << val << " derivative " << dval;
    std::cout << std::endl;
  }
};

const std::string dtcond::condname[numcond] = {"dtMax", "acceleration", "courant", "energy", "smooth", "density", "damage", "sxx", "sxy", "sxz", "syy", "syz"};

class Link {
public:
  int nr, proc;
  void clear() { set(-1, -1); }
  void set(const int &_nr, const int &_proc) { nr = _nr; proc = _proc; }
  template <class T>
  void addParticle(T *part, const int &nr, const int &proc) {
    part->setNext(*this); set(nr, proc);
  }
};

#endif
