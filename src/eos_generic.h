#ifndef SPHLATCH_EOS_GENERIC
#define SPHLATCH_EOS_GENERIC

/*
 *  eos_generic.h
 *
 * 
 *  Created by Andreas Reufer on 26.07.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"
#include "particle_manager.h"
#include "log_manager.h"


namespace sphlatch
{
template<class T_leaftype>
class EOS
{

protected:
typedef sphlatch::ParticleManager partManagerType;
partManagerType& PartManager;
typedef sphlatch::LogManager logManagerType;
logManagerType& Logger;

valvectRefType rho, p, u;
idvectRefType mat;

public:
EOS() :
PartManager(partManagerType::instance()),
Logger(logManagerType::instance()),
rho(PartManager.rho),
p(PartManager.p),
u(PartManager.u),
mat(PartManager.mat)
{
  asLeaf().init();
};

~EOS()
{
};

static EOS<T_leaftype>& instance();
static EOS<T_leaftype>* _instance;

valueType getPressure(const size_t& _i)
{
  return asLeaf().getPressure(_i);
};

valueType getSpeedOfSound(const size_t& _i)
{
  return asLeaf().getSpeedOfSound(_i);
};

void init()
{
  asLeaf().init();
};

private:
T_leaftype& asLeaf()
{
  return static_cast<T_leaftype&>(*this);
}
};

template<class T_leaftype>
EOS<T_leaftype>* EOS<T_leaftype>::_instance = NULL;

template<class T_leaftype>
EOS<T_leaftype>& EOS<T_leaftype>::instance(void)
{
  if (_instance != NULL)
    return *_instance;
  else
    {
      _instance = new EOS;
      return *_instance;
    }
};

}
#endif

