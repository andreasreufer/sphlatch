#include <iostream>
#include <iomanip>

#define SPHLATCH_HDF5
//#undef SPHLATCH_ANEOS_TABLE

#include "typedefs.h"

#include "io_particle.h"
#include "sph_fluid_particle.h"
class particle :
   public sphlatch::SPHfluidPart,
   public sphlatch::energyPart,
   public sphlatch::ANEOSPart,
   public sphlatch::IOPart
{
public:
   ioVarLT getLoadVars()
   {
      ioVarLT vars;

      vars.push_back(storeVar(rho, "rho"));
      vars.push_back(storeVar(mat, "mat"));
      vars.push_back(storeVar(u, "u"));
      vars.push_back(storeVar(p, "p"));
      vars.push_back(storeVar(T, "T"));
      return(vars);
   }

   ioVarLT getSaveVars()
   {
      ioVarLT vars;

      vars.push_back(storeVar(S, "S"));

      return(vars);
   }
};

typedef particle                       partT;

#define SPHLATCH_LOGGER
#include "logger.cpp"
typedef sphlatch::Logger               logT;

#include "eos_super.cpp"
typedef sphlatch::SuperEOS<partT>         eosT;

#include "particle_set.cpp"
typedef sphlatch::ParticleSet<partT>   partSetT;

typedef sphlatch::fType                fType;
typedef sphlatch::iType                iType;

int main(int argc, char* argv[])
{
   eosT& EOS(eosT::instance());
   logT& Logger(logT::instance());

   partSetT parts;

   parts.loadHDF5(argv[1]);

   const size_t nop = parts.getNop();

   Logger.stream << "loaded " << nop << " particles";
   Logger.flushStream();

#ifdef SPHLATCH_ANEOS_TABLE
   //EOS.ANEOS.loadTableU("aneos_tables.hdf5", 4);
   //EOS.ANEOS.loadTableU("aneos_tables.hdf5", 1);
   EOS.ANEOS.loadTableU("aneos_tables.hdf5", 2);
   //EOS.ANEOS.loadTableU("aneos_tables.hdf5", 5);
#endif

   /*size_t bp = 0;
   for (size_t i = 0; i < nop; i++)
   {
      EOS(parts[i]);
      if ( not (parts[i].p == parts[i].p) )
      {
         bp++;
         std::cout << "bad " << i << " " << parts[i].rho << " " << parts[i].u << " " << parts[i].p << " " << parts[i].S << "\n";
      }
   }
   Logger << "pressure";*/

   /*const int mat = 1;
   const fType u = 9.54203e+11;
   const fType rho = 11.9175;*/
   /*const int mat = 2;
   const fType u = 1.67e11    ;
   const fType rho =  3.25572;*/
   
   /*const int mat = 5;
   const fType u = 3.68e11    ;
   const fType rho = 20.25572;*/
   
   const int mat = 5;
   const fType u = 4.72e10    ;
   const fType rho = 13.3841 ;
   fType T, p, S, cs, dummy;
   iType phase;
   
   EOS.aneos.rootU(T, rho, mat, p, u, S, dummy, dummy, dummy, dummy, cs, phase, dummy, dummy);
   std::cout << rho << " " << p << " " << T*11604. << " " << cs << " " << S << " " << phase << "\n";
   
   /*EOS.ANEOS.tableU(T, rho, mat, p, u, S, dummy, dummy, dummy, dummy, cs, phase, dummy, dummy);
   std::cout << rho << " " << p << " " << T*11604. << " " << cs << " " << S << " " << phase << "\n";*/


   /*for (fType rho = 2.0; rho < 15.; rho += 0.1)
   {
   EOS.ANEOS.rootU(T, rho, mat, p, u, S, dummy, dummy, dummy, dummy, cs, phase, dummy, dummy);
   std::cout << rho << " " << p << " " << T*11604. << " " << cs << " " << S << " " << phase << "\n";
   }*/


   //EOS.ANEOS.tableU(T, rho, mat, p, u, S, dummy, dummy, dummy, dummy, cs, phase, dummy, dummy);
   //std::cout << p << " " << T*11604. << " " << cs << " " << S << " " << phase << "\n";
   /*
      partT& p(parts[4]);
      std::cout << p.rho << " " << p.u << " " << p.p << " " << p.T*11700 << "\n";
      EOS(p);
      std::cout << p.T*11700. << " " << p.p << "\n";

   fType T, p, S, dummy;
   iType phase;

   EOS.tableU(T, 2.0, 1, p, 1.e10, S, dummy, dummy, dummy, dummy, dummy, phase,
              dummy,
              dummy);
   Logger << "pressure";*/

   parts.saveHDF5(argv[1]);

   return(EXIT_SUCCESS);
}
