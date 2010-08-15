// some defs

// uncomment for single-precision calculation
//#define SPHLATCH_SINGLEPREC

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include <boost/lexical_cast.hpp>

#include "typedefs.h"
typedef sphlatch::fType fType;
typedef sphlatch::iType iType;

class particle {};

typedef particle partT;

#include "eos_aneos.cpp"
typedef sphlatch::ANEOS<partT> eosT;

int main(int argc, char* argv[])
{
  eosT&      EOS(eosT::instance());

  // in:  rho, u mat
  // out: P, cs, T, phase
  //

  const fType evToK = 11605.333;

  const fType rhoMin = 0.5;
  const fType rhoMax = 4.0;
  const size_t rhoSteps = 100;
  
  const fType Tmin = 1.e2 / evToK;
  const fType Tmax = 1.e4 / evToK;
  const size_t Tsteps = 100;

  const fType rhoLdelta = ( log(rhoMax) - log(rhoMin) ) / rhoSteps;
  const fType TLdelta = ( log(Tmax) - log(Tmin) ) / Tsteps;

  fType P, cs, u;
  
  iType mat, phase;
  std::istringstream matstr(argv[1]);
  matstr >> mat;


  for (size_t i = 0; i < rhoSteps; i++)
  {
    for (size_t j = 0; j < Tsteps; j++)
    {
      const fType rho = rhoMin*exp(rhoLdelta*i);
      const fType T =   Tmin*exp(  TLdelta*j);
  
      std::cout << rho << " " << T*evToK << "\n";
      EOS.getSpecEnergy(rho, T, mat, P, cs, u, phase);

      std::cerr << std::setw(14) << std::setprecision(6) << std::scientific 
                << rho << "\t"
                << T * evToK << "\t"
                << P << "\t"
                << cs << "\t"
                << u << "\t"
                << phase << "\n";
    }
  }

  return EXIT_SUCCESS;
}

