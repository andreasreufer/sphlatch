#ifndef SPH_FLUID_PARTICLE_H
#define SPH_FLUID_PARTICLE_H

/*
 *  sph_fluid_particle.h
 *
 *
 *  Created by Andreas Reufer on 01.03.09
 *  Copyright 2009 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"

namespace sphlatch {
namespace indices
{
enum vect
{
   X, Y, Z
};
};

///
/// basic SPH fluid ghost particle class
///
class SPHfluidGhost {
public:
   fType h, rho, p, cs;
#ifdef SPHLATCH_PADTO64BYTES
private:
   char pad[0];
#endif
};

///
/// basic SPH fluid resident particle class
///
class SPHfluidPart : public SPHfluidGhost
{
  public:
   countsType noneigh;
};

///
/// particle with specific energy
///
class energyGhost { };

class energyPart : public energyGhost
{
  public:
   fType u, dudt;
};

///
/// variable smoothing length ghost
///
class varHGhost { };

class varHPart : public varHGhost
{
  public:
   fType dhdt, divv;
   static cType noneighOpt;
   static fType divvmax;

   void setDivvMax()
   {
     divvmax = divv > divvmax ? divv : divvmax;
   }
};

cType varHPart::noneighOpt;
fType varHPart::divvmax;

///
/// some variables for ANEOS
///
class ANEOSGhost {};

class ANEOSPart : public ANEOSGhost 
{
  public:
  fType T, S;
  iType mat, phase;
};


};
#endif
