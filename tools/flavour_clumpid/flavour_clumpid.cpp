#include <iostream>
#include <sstream>
#include <vector>

//#define SPHLATCH_SINGLEPREC

#include <omp.h>
#define SPHLATCH_OPENMP
#define SPHLATCH_HDF5
#define SPHLATCH_NONEIGH

#include "typedefs.h"
typedef sphlatch::fType     fType;
typedef sphlatch::iType     iType;
typedef sphlatch::cType     cType;
typedef sphlatch::vect3dT   vect3dT;
typedef sphlatch::box3dT    box3dT;

typedef sphlatch::fvectT    fvectT;

const fType finf = sphlatch::fTypeInf;

///
/// define the particle we are using
///
#include "bhtree_particle.h"
#include "clump_particle.h"
#include "io_particle.h"
class particle :
   public sphlatch::treePart,
   public sphlatch::clumpPart,
   /*public sphlatch::SPHfluidPart,
   public sphlatch::energyPart,*/
   public sphlatch::IOPart
{
public:

   ioVarLT getLoadVars()
   {
      ioVarLT vars;

      vars.push_back(storeVar(id, "id"));
      vars.push_back(storeVar(clumpid, "clumpid"));
      return(vars);
   }

   ioVarLT getSaveVars()
   {
      ioVarLT vars;

      vars.push_back(storeVar(clumpidorig, "clumpidorig"));
      
      return(vars);
   }

   iType clumpidorig;
};

typedef particle                               partT;

#include "particle_set.cpp"
typedef sphlatch::ParticleSet<partT>           partSetT;

#include "clump_finder.cpp"
typedef sphlatch::Clumps<partT>                clumpsT;


#include "hdf5_io.cpp"
typedef sphlatch::HDF5File                     H5FT;

// particles are global
partSetT refp, dmpp;
clumpsT  clumps;

class origin {
public:
   iType clumpid;
};

typedef origin  originT;
typedef origin& originRT;

int main(int argc, char* argv[])
{
   if (not (argc == 3) )
   {
      std::cerr <<
      "usage: flavour_clumpid <refdump> <file>\n";
      return(1);
   }

   std::string refname = argv[1];
   std::string dmpname = argv[2];

   // load the reference particles
   refp.loadHDF5(refname);
   const size_t norp = refp.getNop();

   std::map<iType, origin> originMap;

   for (size_t i = 0; i < norp; i++)
     originMap[ refp[i].id ].clumpid = refp[i].clumpid;

   
   dmpp.loadHDF5(dmpname);
   const size_t nodp = dmpp.getNop();
   for (size_t i = 0; i < nodp; i++)
     dmpp[i].clumpidorig = originMap[ dmpp[i].id ].clumpid;
   
   dmpp.attributes["reftime"] = refp.attributes["time"];
   dmpp.saveHDF5(dmpname);

   return(0);
}
