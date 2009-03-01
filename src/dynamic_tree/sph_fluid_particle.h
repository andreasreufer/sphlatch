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
///
/// basic SPH fluid ghost particle class
///
class SPHfluidGhost {
public:
   fType        h, rho, p, cs;
#ifdef SPHLATCH_PADTO64BYTES
private:
   char pad[0];
#endif
};

///
/// basic SPH fluid resident particle class
///
class SPHfluidPart : public SPHfluidGhost {
  countsType noneigh;
};

///
/// particle with specific energy
///
class energyGhost {};

class energyPart : public energyGhost
{
  fType u, dudt;
};

///
/// variable smoothing length ghost
///
class varHGhost {};

class varHPart : public varHGhost
{
  fType dhdt, divv;
};

};
#endif
