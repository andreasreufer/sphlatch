// some defs

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#define SPHLATCH_HDF5

#include "typedefs.h"
typedef sphlatch::fType   fType;
typedef sphlatch::iType   iType;

#include "bhtree_particle.h"
#include "sph_fluid_particle.h"
#include "io_particle.h"
class particle :
   public sphlatch::treePart,
   public sphlatch::SPHfluidPart,
   public sphlatch::IOPart
{
public:
   ioVarLT getLoadVars()
   {
      ioVarLT vars;

      return(vars);
   }

   ioVarLT getSaveVars()
   {
      ioVarLT vars;

      vars.push_back(storeVar(pos, "pos"));
      vars.push_back(storeVar(id, "id"));
      vars.push_back(storeVar(m, "m"));
      vars.push_back(storeVar(h, "h"));
      return(vars);
   }
};
typedef particle                       partT;

#include "particle_set.cpp"
typedef sphlatch::ParticleSet<partT>   partSetT;

#include "lattice_hcp.h"
typedef sphlatch::LatticeHCP           lattT;

int main(int argc, char* argv[])
{
   if (not (argc == 3))
   {
      std::cerr <<
      "usage: vanilla_sphere <dump> <nop>\n";
      return(1);
   }
   partSetT parts;

   std::string outFile = argv[1];

   std::stringstream nopstr(argv[2]);
   size_t            nopw;
   nopstr >> nopw;

   ///
   /// 1.03 is a fudge factor :-)
   ///
   const fType rMax          = 1.;
   const fType fillingFactor = 0.740;
   const fType latticeLength = 2. * rMax / pow(fillingFactor * nopw
                                               / 1.03, 1. / 3.);

   ///
   /// now place the SPH particles on a lattice
   ///
   lattT Lattice(latticeLength, 1.1 * rMax, 1.1 * rMax, 1.1 * rMax);

   size_t pc = 0;
   Lattice.first();
   while (!Lattice.isLast)
   {
      if (Lattice.rCur < rMax)
         pc++;
      Lattice.next();
   }
   std::cerr << "you asked for " << nopw
             << " and you get " << pc << " particles\n";

   parts.resize(pc);

   const fType totV      = (4. * M_PI / 3.) * pow(rMax, 3.);
   const fType partV     = totV / pc;
   const fType smoLength = 0.85 * latticeLength;

   std::cerr << " particle volume: " << partV << "\n"
             << " total    volume: " << totV << "\n";

   pc = 0;
   Lattice.first();
   while (!Lattice.isLast)
   {
      if (Lattice.rCur < rMax)
      {
         parts[pc].id     = pc;
         parts[pc].pos[0] = Lattice.xCur;
         parts[pc].pos[1] = Lattice.yCur;
         parts[pc].pos[2] = Lattice.zCur;

         ///
         /// the particle volume is stored in the mass variable
         ///
         /// so, multiplied with the density, one gets the
         /// particle mass
         ///
         parts[pc].h = smoLength;
         parts[pc].m = partV;

         pc++;
      }
      Lattice.next();
   }

   std::cerr << " -> " << outFile << "\n";
   parts.step = 0;
   parts.saveHDF5(outFile);

   return(EXIT_SUCCESS);
}
