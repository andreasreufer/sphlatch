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
      std::cerr << "rootPX <matid> <T [K]> <p [barye]>\n";
      return(EXIT_FAILURE);
   }

   eosT& EOS(eosT::instance());

   std::istringstream matStr(argv[1]);
   iType mat;
   matStr >> mat;

   std::istringstream   TStr(argv[2]);
   fType TK;
     TStr >> TK;

   std::istringstream   pStr(argv[3]);
   fType p ;
     pStr >> p;

   const fType eVinK = 11604.505;

   fType rho, u, S, cv, dummy, cs;
   iType phase;

   EOS.rootP(TK / eVinK, rho, mat, p, u, S, cv, dummy, dummy, dummy, cs, phase,
             dummy,
             dummy);
   
   std::cout << "rho:   " << rho << "\n" 
             << "u:     " << u << "\n"
             << "S:     " << S << "\n"
             << "cv:    " << cv << "\n"
             << "cs:    " << cs << "\n"
             << "phase: " << phase << "\n";
   return(EXIT_SUCCESS);
}
