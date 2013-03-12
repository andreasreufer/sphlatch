// some defs

// uncomment for single-precision calculation
//#define SPHLATCH_SINGLEPREC

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

//#include <boost/lexical_cast.hpp>

#define SPHLATCH_HDF5

#include "typedefs.h"
typedef sphlatch::fType   fType;

const fType finf = sphlatch::fTypeInf;

#include "io_particle.h"

class particle : public sphlatch::IOPart {
public:
   fType Tmax, uav, u;

   ioVarLT getLoadVars()
   {
      ioVarLT vars;

      vars.push_back(storeVar(u, "u"));

      return(vars);
   }

   ioVarLT getSaveVars()
   {
      ioVarLT vars;

      vars.push_back(storeVar(Tmax, "Tmax"));
      vars.push_back(storeVar(uav, "uav"));

      return(vars);
   }
};

typedef particle                       partT;

#include "particle_set.cpp"
typedef sphlatch::ParticleSet<partT>   partSetT;


int main(int argc, char* argv[])
{
   std::string fname(argv[1]);

   partSetT parts;

   parts.loadHDF5(fname);
   const size_t nop = parts.getNop();

   for (size_t i = 0; i < nop; i++)
   {
     parts[i].uav = 0.;
     parts[i].Tmax = -finf;
   }
   parts.attributes["ethermadded"] = 0.;
   parts.saveHDF5(fname);

   return(EXIT_SUCCESS);
}
