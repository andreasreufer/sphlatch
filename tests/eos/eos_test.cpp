#include <iostream>
#include <iomanip>

#define SPHLATCH_HDF5

#include "typedefs.h"

#include "sph_fluid_particle.h"
class particle :
  public sphlatch::SPHfluidPart,
  public sphlatch::energyPart,
  public sphlatch::ANEOSPart
{
  public:
    int id;
};

typedef particle partT;

#include "eos_aneos.cpp"
typedef sphlatch::ANEOS<partT> eosT;

typedef sphlatch::fType fType;
typedef sphlatch::iType iType;

int main(int argc, char* argv[])
{
  eosT& EOS(eosT::instance());  


#ifdef SPHLATCH_ANEOS_TABLE  
  //EOS.storeTableU("aneos_tables.hdf5",1);
#endif
  
  //EOS.rootu(4.0, 1.e10, 4);

  //EOS.rootS(4.0, 3.78352e10, 4);
  /*EOS.rootS(3.20, 3.00000e11, 4);
  EOS.rootS(3.32, 3.00000e11, 4);
  EOS.rootS(3.50, 3.00000e11, 4);
  EOS.rootS(4.00, 3.00000e11, 4);
  EOS.rootS(4.50, 3.00000e11, 4);*/
  fType T, rho, p, u, S, cv, dpdt, dpdr, fkros, cs, fme, fma;
  iType mat, kpa;
  
  const fType eVinK = 11605.;
  const fType erggeVinJkgK = 1.1605e8;

#ifdef SPHLATCH_ANEOS
  mat = 2;
#else
  mat = 1;
#endif


  partT part;
  part.rho = 4.0;
  part.u   = 1.e10;
  part.mat = mat;
  
  


  const fType pref = 1.e6;

  for(fType TK = 200.; TK < 50.e3; TK +=  30.)
  //for(fType TK = 100.; TK < 50.e2; TK +=  10.)
  {
    EOS.rootP(TK/eVinK, rho, mat, pref, u, S, cv, dpdt, dpdr, fkros, cs, kpa, fme, fma);
    //std::cout << TK << " " << rho << " " << u << " " << S/erggeVinJkgK << " " << cs << "\n";
  }

  std::cout << "          rho     u           S            p       cs     T\n";

  fType TK = 1000.;
  // T,p -> rho, u, S
  EOS.rootP(TK/eVinK, rho, mat, pref, u, S, cv, dpdt, dpdr, fkros, cs, kpa, fme, fma);
  std::cout << "rootP:    " << rho << " " << u << " " << S << " " << pref << " " << cs << " " << TK << "\n";

  const fType Sref = S;
  const fType uref = u;

  // rho, S -> p, u
  EOS.rootS(T, rho, mat, p, u, S, cv, dpdt, dpdr, fkros, cs, kpa, fme, fma);
  std::cout << "rootS:    " << rho << " " << u << " " << S << " " << p << " " << cs << " " << T*eVinK << "\n";

#ifdef SPHLATCH_ANEOS_TABLE  
  EOS.loadTableS("aneos_tables.hdf5",mat);
  std::cout << "loaded table!\n";
#endif
  
  u = uref;
  
  // rho, S -> p, u
  EOS.tableS(T, rho, mat, p, u, S, cv, dpdt, dpdr, fkros, cs, kpa, fme, fma);
  std::cout << "tableS:   " << rho << " " << u << " " << S << " " << p << " " << cs << " " << T*eVinK << "\n";

#ifdef SPHLATCH_ANEOS_TABLE  
  //EOS.storeTableS("aneos_tables.hdf5",mat);
#endif

  // rho, u -> p, S
  EOS.rootU(T, rho, mat, p, u, S, cv, dpdt, dpdr, fkros, cs, kpa, fme, fma);
  std::cout << "rootU:    " << rho << " " << u << " " << S << " " << p << " " << cs << " " << T*eVinK << "\n";

#ifdef SPHLATCH_ANEOS_TABLE  
  EOS.loadTableU("aneos_tables.hdf5",mat);
#endif
  // rho, u -> p, S
  EOS.tableU(T, rho, mat, p, u, S, cv, dpdt, dpdr, fkros, cs, kpa, fme, fma);
  std::cout << "tableU:   " << rho << " " << u << " " << S << " " << p << " " << cs << " " << T*eVinK << "\n";

  part.rho = rho;
  part.u   = u;
  EOS(part);
  std::cout << part.S << " " << part.p << " " << part.cs << " " << part.T*11700. << " " << part.phase << "\n";
  
  part.rho = 1.e-11;
  EOS(part);
  std::cout << part.S << " " << part.p << " " << part.cs << " " << part.T*11700. << " " << part.phase << "\n";

  part.rho = 100.;
  part.u = 1.e15;
  EOS(part);
  std::cout << part.S << " " << part.p << " " << part.cs << " " << part.T*11700. << " " << part.phase << "\n";

#ifdef SPHLATCH_ANEOS_TABLE  
  //EOS.storeTableU("aneos_tables.hdf5",mat);
#endif
  
  return EXIT_SUCCESS;
}
