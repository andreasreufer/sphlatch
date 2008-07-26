#ifndef SPHLATCH_EOS_TILLOTSON
#define SPHLATCH_EOS_TILLOTSON

/*
 *  eos_tillotson.h
 *
 *
 *  Created by Andreas Reufer on 26.07.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"
#include "eos_generic.h"

namespace sphlatch {

class Tillotson : public EOS<Tillotson> {

public:

void init()
{
  Logger << "init Tillotson EOS";
};

valueType getPressure(const size_t& _i)
{
/*  const valueType gamma = PartManager.attributes["gamma"];

  valvectRefType p( PartManager.p);
  valvectRefType rho( PartManager.rho);
  valvectRefType u( PartManager.u);
  
  p = (gamma - 1) * (boost::numeric::ublas::element_prod(u, rho));*/
  return 0.;
};

valueType getSpeedOfSound(const size_t& _i)
{
  return 0.;
};

valueType getTemperature(const size_t& _i)
{
  return 0.;
};

private:

};
}
#endif
