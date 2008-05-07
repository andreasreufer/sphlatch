#ifndef SPHLATCH_PARTICLE_MANAGER_H
#define SPHLATCH_PARTICLE_MANAGER_H

#include <map>
#include <set>
#include <boost/assign/std/set.hpp>
/* #include <boost/assign/std/vector.hpp>*/

#include "typedefs.h"

namespace sphlatch
{

enum vectorialIndices { X, Y, Z };

class ParticleManager
{
public:
typedef ParticleManager self_type;
typedef ParticleManager& self_reference;
typedef ParticleManager* self_pointer;

typedef std::map< std::string, valvectPtrType > valVectPtrMap;
typedef std::map< std::string, matrixPtrType > matrixPtrMap;
typedef std::map< std::string, idvectPtrType > idVectPtrMap;

typedef std::set< std::string > stringSet;

static self_reference instance(void);

///
/// methods to set up the vars needed
///
void useBasicSPH(void);
void useEnergy(void);
void useTimedepEnergy(void);
void useEntropy(void);
void useTimedepEntropy(void);
void useTimedepH(void);
void useGravity(void);

void useAVMonaghan(void);
void useAVBalsara(void);
void useAVHernquist(void);
void useAVTimedepAlpha(void);

void setNoParts(size_t _noParts, size_t _noGhostParts, bool _keepData);
void setNoParts(size_t _noParts, size_t _noGhostParts);
void setNoParts(size_t _noParts);

///
/// get the name of a int/scalar/vect value reference
///
std::string getName(matrixRefType _matrixRef);
std::string getName(valvectRefType _valvectRef);
std::string getName(idvectRefType _idvectRef);

///
/// get the int/scalar/vect value reference for a name
///
matrixPtrType  getVectRef(std::string _name);
valvectPtrType getScalarRef(std::string _name);
idvectPtrType  getIdRef(std::string _name);

///
/// vector quantities
///
matrixType pos, vel, acc, rotv, M, I, S;

///
/// scalar quantities
///
valvectType m, h, dhdt, rho, drhodt, p, u, dudt, A, dAdt, alpha, dalphadt,
            divv, mumax, q, eps;

///
/// integers
///
idvectType id, noneigh;

///
/// the integration step
///
size_t step;

protected:
ParticleManager(void);
~ParticleManager(void);

private:
///
/// the quantities used by local particles
///
valVectPtrMap usedScalars;
matrixPtrMap usedVectors;
idVectPtrMap usedIntegers;

///
/// number of local/ghost particles
///
private:
size_t noLocalParts, noGhostParts;

public:
size_t getNoLocalParts();
size_t getNoGhostParts();

///
/// the quantities which are ghosts, when used
///
private:
stringSet isGhostVarSet;
bool isGhostVar(std::string _searchString);

///
/// attribute stuff
///
public:
attrMapType attributes;
bool attrExists(std::string _name);

static self_pointer _instance;
};

ParticleManager::self_pointer ParticleManager::_instance = NULL;

ParticleManager::self_reference ParticleManager::instance(void)
{
  if (_instance != NULL)
    return *_instance;
  else
    {
      _instance = new ParticleManager;
      return *_instance;
    }
}

ParticleManager::ParticleManager(void)
{
  ///
  /// define the amount of variables needed in every case
  /// (position and ID)
  ///
  usedVectors[ "pos" ] = &pos;
  pos.resize(0,3);

  usedIntegers[ "id" ] = &id;

  using namespace boost::assign;
  ///
  /// define the quantities used for ghosts
  /// vectorial quantities:
  ///
  isGhostVarSet += "pos", "vel", "rotv";
  /// scalar quantities:
  isGhostVarSet += "m", "h", "rho", "p", "u", "A", "alpha", "divv";
  /// integers:
  isGhostVarSet += "id";

  step = 0;
}


ParticleManager::~ParticleManager(void)
{
}

///
/// everything for nbody
///
void ParticleManager::useGravity()
{
  usedVectors[ "vel" ] = &vel;
  usedVectors[ "acc" ] = &acc;
  
  // 3D
  vel.resize(vel.size1(), 3);
  acc.resize(acc.size1(), 3);

  usedScalars[ "m" ] = &m;
  usedScalars[ "eps" ] = &eps;
}

///
/// most basic SPH variables
///
void ParticleManager::useBasicSPH()
{
  usedVectors[ "vel" ] = &vel;
  usedVectors[ "acc" ] = &acc;

  // 3D
  vel.resize(vel.size1(), 3);
  acc.resize(acc.size1(), 3);
  
  usedScalars[ "m" ] = &m;
  usedScalars[ "h" ] = &h;
  usedScalars[ "rho" ] = &rho;
  usedScalars[ "p" ] = &p;
}

///
/// use internal energy
///
void ParticleManager::useEnergy()
{
  usedScalars[ "u" ] = &u;
}

///
/// time dependent internal energy
///
void ParticleManager::useTimedepEnergy()
{
  usedScalars[ "dudt" ] = &dudt;
}

///
/// use entropy
///
void ParticleManager::useEntropy()
{
  usedScalars[ "A" ] = &A;
}

///
/// time dependent entropy
///
void ParticleManager::useTimedepEntropy()
{
  usedScalars[ "dAdt" ] = &dAdt;
}

///
/// time dependent smoothing length
///
void ParticleManager::useTimedepH()
{
  usedScalars[ "dhdt" ] = &dhdt;
  usedIntegers[ "noneigh" ] = &id;
}

///
/// vars for Monaghan (1989) artificial viscosity
///
void ParticleManager::useAVMonaghan()
{
  usedScalars[ "mumax" ] = &mumax;
}

///
/// time dependent alpha
///
void ParticleManager::useAVTimedepAlpha()
{
  usedScalars[ "alpha" ] = &alpha;
  usedScalars[ "dalphadt" ] = &dalphadt;
}

///
/// balsara switch
///
void ParticleManager::useAVBalsara()
{
  usedVectors[ "rotv" ] = &rotv;
  
  // 3D
  rotv.resize(rotv.size1(), 3);
}

///
/// vars for Hernquist&Katz artificial viscosity
///
void ParticleManager::useAVHernquist()
{
  usedScalars[ "q" ] = &q;
}

void ParticleManager::setNoParts(size_t _noLocalParts, size_t _noGhostParts, 
                                 bool _keepData)
{
  noLocalParts = _noLocalParts;
  noGhostParts = _noGhostParts;

  size_t locSize = noLocalParts;
  size_t ghostSize = noLocalParts + noGhostParts;

  matrixPtrMap::iterator vectorsItr = usedVectors.begin();

  while (vectorsItr != usedVectors.end())
    {
      if (isGhostVar(vectorsItr->first))
        {
          vectorsItr->second->resize(ghostSize,
                                     vectorsItr->second->size2(), _keepData);
        }
      else
        {
          vectorsItr->second->resize(locSize,
                                     vectorsItr->second->size2(), _keepData);
        }
      vectorsItr++;
    }

  valVectPtrMap::iterator scalarsItr = usedScalars.begin();
  while (scalarsItr != usedScalars.end())
    {
      if (isGhostVar(scalarsItr->first))
        {
          scalarsItr->second->resize(ghostSize, _keepData);
        }
      else
        {
          scalarsItr->second->resize(locSize, _keepData);
        }
      scalarsItr++;
    }

  idVectPtrMap::iterator intsItr = usedIntegers.begin();
  while (intsItr != usedIntegers.end())
    {
      if (isGhostVar(intsItr->first))
        {
          intsItr->second->resize(ghostSize, _keepData);
        }
      else
        {
          intsItr->second->resize(locSize, _keepData);
        }
      intsItr++;
    }
}

///
/// specialization in case you want to keep the data in the containers
///
void ParticleManager::setNoParts(size_t _noLocalParts, size_t _noGhostParts)
{
  setNoParts(_noLocalParts, _noGhostParts, true);
}

///
/// specialization in case you don't need ghosts and data can be discarded
///
void ParticleManager::setNoParts(size_t _noLocParts)
{
  setNoParts(_noLocParts, 0, false);
}

size_t ParticleManager::getNoLocalParts()
{
  return noLocalParts;
}

size_t ParticleManager::getNoGhostParts()
{
  return noGhostParts;
}

///
/// is the variable with a given name used by ghosts?
///
bool ParticleManager::isGhostVar(std::string _searchString)
{
  stringSet::const_iterator strItr = isGhostVarSet.begin();

  strItr = isGhostVarSet.find(_searchString);
  return(strItr != isGhostVarSet.end());
}

///
/// get a vectorial quantitiy references name
///
std::string ParticleManager::getName(matrixRefType _matrixRef)
{
  matrixPtrMap::iterator vectorsItr = usedVectors.begin();

  while (vectorsItr != usedVectors.end())
    {
      if (vectorsItr->second == &_matrixRef)
        {
          return vectorsItr->first;
        }
      vectorsItr++;
    }
  return "";
}

///
/// get a scalar quantitiy references name
///
std::string ParticleManager::getName(valvectRefType _valvectRef)
{
  valVectPtrMap::iterator scalarsItr = usedScalars.begin();

  while (scalarsItr != usedScalars.end())
    {
      if (scalarsItr->second == &_valvectRef)
        {
          return scalarsItr->first;
        }
      scalarsItr++;
    }
  return "";
}

///
/// get an integer quantitiy references name
///
std::string ParticleManager::getName(idvectRefType _idvectRef)
{
  idVectPtrMap::iterator intsItr = usedIntegers.begin();

  while (intsItr != usedIntegers.end())
    {
      if (intsItr->second == &_idvectRef)
        {
          return intsItr->first;
        }
      intsItr++;
    }
  return "";
}

///
/// get a vectorial quantity pointer by name
///
matrixPtrType ParticleManager::getVectRef(std::string _name)
{
  matrixPtrMap::iterator matrItr = usedVectors.find( _name );

  if ( matrItr == usedVectors.end() )
  {
    return NULL;
  }
  else
  {
    return matrItr->second;
  }
}

///
/// get a scalar quantity pointer by name
///
valvectPtrType ParticleManager::getScalarRef(std::string _name)
{
  valVectPtrMap::iterator vectItr = usedScalars.find( _name );

  if ( vectItr == usedScalars.end() )
  {
    return NULL;
  }
  else
  {
    return vectItr->second;
  }
}

///
/// get an integer quantity pointer by name
///
idvectPtrType ParticleManager::getIdRef(std::string _name)
{
  idVectPtrMap::iterator intItr = usedIntegers.find( _name );

  if ( intItr == usedIntegers.end() )
  {
    return NULL;
  }
  else
  {
    return intItr->second;
  }
}

///
/// does the attribute exist?
///
bool ParticleManager::attrExists(std::string _name)
{
  return ( attributes.find( _name ) != attributes.end() );
}

};

#endif
