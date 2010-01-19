#ifndef SPHLATCH_CONSTANTS
#define SPHLATCH_CONSTANTS

/*
 *  constants.h
 *
 *
 *  Created by Andreas Reufer on 21.12.09
 *  Copyright 2009 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"

namespace sphlatch {
namespace constants {
// eV expressed in Kelvin
const fType eVinK = 11604.505;

namespace unitsSI
{
const fType G = 6.67428e-11;   // m^3 kg^-1 s^-2

const fType M_earth = 5.9736e24;   // kg
const fType R_earth = 6.3781e6;    // m
};

namespace unitsCGS
{
// conversion constants to derive quantities from SI
const fType unitLength = 1.e2; // cm in m
const fType unitTime   = 1.e0; // s  in s
const fType unitMass   = 1.e3; // g  in kg

const fType G = unitsSI::G * (unitLength * unitLength * unitLength)
                / (unitMass * unitTime * unitTime);

const fType M_earth = unitsSI::M_earth * unitMass;
const fType R_earth = unitsSI::R_earth * unitLength;
};
};
};

#endif
