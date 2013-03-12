#include <iostream>
#include <iomanip>
#include "typedefs.h"

class particle { };

#include "eos_aneos.cpp"
typedef sphlatch::ANEOS<particle>   eosT;

typedef sphlatch::fType             fType;
typedef sphlatch::iType             iType;

int main(int argc, char* argv[])
{
   if (argc != 4)
   {
      std::cerr << "rootSX <matid> <rho [g/cm^3]> <S [erg/g/K]>\n";
      return(EXIT_FAILURE);
   }

   eosT& EOS(eosT::instance());

   std::istringstream matStr(argv[1]);
   iType mat;
   matStr >> mat;

   std::istringstream   rhoStr(argv[2]);
   fType rho;
     rhoStr >> rho;

   std::istringstream   SStr(argv[3]);
   fType S;
     SStr >> S;

   const fType eVinK = 11604.505;
   S *= eVinK;

   fType p, T, u, cv, dummy, cs, rhoL, rhoH;
   iType phase;

   EOS.rootS(T, rho, mat, p, u, S, cv, dummy, dummy, dummy, cs, phase,
             rhoL, rhoH);
   
   std::cout << "p [barye]:     " << p   << "\n" 
             << "T [K]:         " << T*eVinK << "\n"
             //<< "S: [erg/g/eV]  " << S << "\n"
             << "u: [erg/g]    " << u << "\n"
             << "cv: [erg/g/ev] " << cv << "\n"
             << "cs: [cm/s]     " << cs << "\n"
             << "phase:         " << phase << "\n"
             << "rho low:       " << rhoL << "\n"
             << "rho high:      " << rhoH << "\n";

   return(EXIT_SUCCESS);
}
