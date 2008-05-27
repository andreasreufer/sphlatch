#ifndef BHTREE_INTEGRATOR_GENERIC_H
#define BHTREE_INTEGRATOR_GENERIC_H

/*
 *  integrator_generic.h
 *
 *
 *  Created by Andreas Reufer on 26.05.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

//#include "communication_manager.h"
#include "typedefs.h"
#include "particle_manager.h"

namespace sphlatch {

template<class T_leaftype>

class MetaIntegrator {

//typedef sphlatch::CommunicationManager commManagerType;
typedef sphlatch::ParticleManager      partManagerType;

private:
T_leaftype& asLeaf()
{
  return static_cast<T_leaftype&>(*this);
}

public:
///
/// \brief constructor:
/// 
///
MetaIntegrator() :
//CommManager(commManagerType::instance()),
PartManager(partManagerType::instance())
{
};

///
/// destructor:
///
~MetaIntegrator(void)
{
  /// do not forget to unregister the integrating vars!
}

protected:
//commManagerType& CommManager;
partManagerType& PartManager;

///
/// variables
///

///
/// register variables for 1st order integration
///
public:
void regIntegration( matrixRefType _var,
                     matrixRefType _devVar )
{
  asLeaf().regIntegration(_var, _devVar );
}

void regIntegration( valvectRefType _var,
                     valvectRefType _devVar )
{
  asLeaf().regIntegration(_var, _devVar );
}

///
/// register variables for 2nd order integration
///
void regIntegration( matrixRefType _var,
                     matrixRefType _devVar,
                     matrixRefType _ddevVar )
{
  asLeaf().regIntegration(_var, _devVar, _ddevVar );
}

void regIntegration( valvectRefType _var,
                     valvectRefType _devVar,
                     valvectRefType _ddevVar )
{
  asLeaf().regIntegration(_var, _devVar, _ddevVar );
}

///


};

};

#endif
