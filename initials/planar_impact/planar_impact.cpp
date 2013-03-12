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
   public sphlatch::movingPart,
   public sphlatch::SPHfluidPart,
   public sphlatch::IOPart,
   public sphlatch::ANEOSPart,
   public sphlatch::energyPart
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
      vars.push_back(storeVar(vel, "vel"));
      vars.push_back(storeVar(id, "id"));
      vars.push_back(storeVar(m, "m"));
      vars.push_back(storeVar(h, "h"));
      vars.push_back(storeVar(u, "u"));
      vars.push_back(storeVar(rho, "rho"));
      vars.push_back(storeVar(mat, "mat"));
      return(vars);
   }
};
typedef particle                       partT;

#include "eos_super.cpp"
typedef sphlatch::SuperEOS<partT>      eosT;

#include "particle_set.cpp"
typedef sphlatch::ParticleSet<partT>   partSetT;

#include "lattice_hcp.h"
typedef sphlatch::LatticeHCP           lattT;

int main(int argc, char* argv[])
{
   if (not (argc == 9))
   {
      std::cerr <<
      "usage: planar_impact <dump> <nop> <mat1> <T1> <mat2> <T2> <impvel> <boxlength>\n";
      return(1);
   }
   partSetT parts;

   std::string outFile = argv[1];

   std::stringstream nopstr(argv[2]);
   size_t            nopw;
   nopstr >> nopw;

   std::stringstream mat1str(argv[3]);
   iType             mat1;
   mat1str >> mat1;

   std::stringstream T1str(argv[4]);
   fType             T1;
   T1str >> T1;

   std::stringstream mat2str(argv[5]);
   iType             mat2;
   mat2str >> mat2;

   std::stringstream T2str(argv[6]);
   fType             T2;
   T2str >> T2;

   std::stringstream impvelstr(argv[7]);
   fType             impvel;
   impvelstr >> impvel;

   std::stringstream lstr(argv[8]);
   fType             l;
   lstr >> l;


   ///
   /// 1.05 is a fudge factor :-)
   ///
   // volume and lattice length for one of the bodies
   const fType fillingFactor = 0.740;
   const fType totV = 4.*l*l*l;
   const fType latl = 2.*pow( fillingFactor * totV / (nopw) / 1.05, 1./3.);
   //const fType smol = 0.85 * latl;
   const fType smol = 1.15 * latl;

   eosT& EOS(eosT::instance());

   const fType eVinK = 11604.505;
   fType       rho1, rho2, u1, u2, dummy;
   iType       idummy;
   EOS.aneos.rootP(T1 / eVinK, rho1, mat1, 1.e3, u1, dummy, dummy, dummy, dummy,
                   dummy, dummy, idummy, dummy, dummy);
   EOS.aneos.rootP(T2 / eVinK, rho2, mat2, 1.e3, u2, dummy, dummy, dummy, dummy,
                   dummy, dummy, idummy, dummy, dummy);

   std::cout << "mat1 (target):   rho " << rho1 << "  u " << u1 << "\n";
   std::cout << "mat2 (impactor): rho " << rho2 << "  u " << u2 << "\n";


   ///
   /// now place the SPH particles on a lattice
   ///
   lattT Lattice(latl, 1.01 * l, 1.01 * l, 0.51 * l);


   size_t pc = 0;
   Lattice.first();
   while (!Lattice.isLast)
   {
      pc++;
      Lattice.next();
   }

   pc *= 2;
   std::cerr << "you asked for " << nopw
             << " and you get " << pc << " particles\n";

   parts.resize(pc);

   const fType pV        = 2.*totV / pc;

   const fType z0dist = 5. * smol;

   std::cerr << " particle volume: " << pV << "\n"
             << " total    volume: " << 2.*totV << "\n";

   const size_t po = pc / 2;
   pc = 0;

   Lattice.first();
   while (!Lattice.isLast)
   {
      // the resting target plate
      parts[pc].id     = pc;
      parts[pc].pos[0] = Lattice.xCur;
      parts[pc].pos[1] = Lattice.yCur;
      parts[pc].pos[2] = Lattice.zCur - 0.50 * l;

      parts[pc].vel[0] = 0.;
      parts[pc].vel[1] = 0.;
      parts[pc].vel[2] = 0.;

      parts[pc].h   = smol;
      parts[pc].m   = pV * rho1;
      parts[pc].rho = rho1;
      parts[pc].u   = u1;

      // the impacting plate
      parts[pc + po].id     = pc + 2000000;
      parts[pc + po].pos[0] = Lattice.xCur;
      parts[pc + po].pos[1] = Lattice.yCur;
      parts[pc + po].pos[2] = Lattice.zCur + z0dist + 0.50 * l;

      parts[pc + po].vel[0] = 0.;
      parts[pc + po].vel[1] = 0.;
      parts[pc + po].vel[2] = -impvel;

      parts[pc + po].h   = smol;
      parts[pc + po].m   = pV * rho2;
      parts[pc + po].rho = rho2;
      parts[pc + po].u   = u2;


      ///
      /// the particle volume is stored in the mass variable
      ///
      /// so, multiplied with the density, one gets the
      /// particle mass
      ///
      //parts[pc].h = smol;
      //parts[pc].m = pV;

      pc++;
      Lattice.next();
   }

   // set velocity

   std::cerr << " -> " << outFile << "\n";
   parts.step = 0;
   parts.attributes["time"] = -(z0dist - latl) / impvel;
   parts.saveHDF5(outFile);

   return(EXIT_SUCCESS);
}
