#ifndef SPHLATCH_PARTICLE_MANAGER_H
#define SPHLATCH_PARTICLE_MANAGER_H

#include <map>
#include <set>
#include <boost/assign/std/set.hpp>

#include "typedefs.h"

namespace sphlatch
{

namespace vectindices
{
enum Vectindices { X, Y, Z };
}

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
void useIntegratedRho(void);

void useAVMonaghan(void);
void useAVBalsara(void);
void useAVHernquist(void);
void useAVTimedepAlpha(void);

void useMaterials(void);
void usePhase(void);

void setNoParts(size_t _noParts, size_t _noGhostParts);
void setNoParts(size_t _noParts);

void addParts(size_t _addNoParts);

void resizeAll(bool _keep);
void resizeAll(void);

void resize(matrixRefType _matrixRef, bool _keep);
void resize(valvectRefType _valvectRef, bool _keep);
void resize(idvectRefType _idvectRef, bool _keep);

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
/// register an auxilliary quantity, for example integration vars
///
void regQuantity(matrixRefType _matrixRef,   std::string _name);
void regQuantity(valvectRefType _valvectRef, std::string _name);
void regQuantity(idvectRefType _idvectRef,   std::string _name);

///
/// unregister a quantity again, for example when an integrator is destructed
///
void unRegQuantity(matrixRefType  _matRef);
void unRegQuantity(valvectRefType _valvectRef);
void unRegQuantity(idvectRefType  _idvectRef);

///
/// vector quantities
///
matrixType pos, vel, acc, rotv, M, I, S;

///
/// scalar quantities
///
valvectType m, h, dhdt, rho, drhodt, p, u, dudt, A, dAdt, alpha, dalphadt,
            divv, mumax, q, eps, dtav, dt, cs;

///
/// integers
///
idvectType id, noneigh, mat, phase;

///
/// blacklist for particles to be deleted
///
bitsetType blacklisted;

///
/// the integration step and the substep of the integrator
///
int step, substep;

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
size_t getNoTotalParts();

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
  pos.resize(0, 3);

  usedIntegers[ "id" ] = &id;

  using namespace boost::assign;
  ///
  /// define the quantities used for ghosts
  /// vectorial quantities:
  ///
  isGhostVarSet += "pos", "vel", "rotv";
  /// scalar quantities:
  isGhostVarSet += "m", "h", "rho", "p", "A",
                    "alpha", "divv", "eps", "cs";
  /// integers:
  isGhostVarSet += "id";

  step = 0, substep = 0;
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
  usedScalars[ "cs" ] = &cs;
  
  usedIntegers[ "noneigh" ] = &noneigh;
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
  usedScalars[ "divv" ] = &divv;
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

///
/// the material variable
///
void ParticleManager::useMaterials()
{
  usedIntegers[ "mat" ] = &mat;
}

///
/// the material variable
///
void ParticleManager::usePhase()
{
  usedIntegers[ "phase" ] = &phase;
}

///
/// the material variable
///
void ParticleManager::useIntegratedRho()
{
  usedScalars[ "drhodt" ] = &drhodt;
}

void ParticleManager::setNoParts(size_t _noLocalParts, size_t _noGhostParts)
{
  noLocalParts = _noLocalParts;
  noGhostParts = _noGhostParts;
}

///
/// specialization in case you don't need ghosts
///
void ParticleManager::setNoParts(size_t _noLocParts)
{
  setNoParts(_noLocParts, 0);
}

///
///
//
void ParticleManager::addParts(size_t _addNoParts)
{
  setNoParts(getNoLocalParts() + _addNoParts, 0);
  resizeAll(true);
}

///
/// resize ALL quantities
///
void ParticleManager::resizeAll(bool _keep)
{
  matrixPtrMap::iterator vectorsItr = usedVectors.begin();

  while (vectorsItr != usedVectors.end())
    {
      resize( *(vectorsItr->second), _keep);
      vectorsItr++;
    }

  valVectPtrMap::iterator scalarsItr = usedScalars.begin();
  while (scalarsItr != usedScalars.end())
    {
      resize( *(scalarsItr->second), _keep);
      scalarsItr++;
    }

  idVectPtrMap::iterator intsItr = usedIntegers.begin();
  while (intsItr != usedIntegers.end())
    {
      resize( *(intsItr->second), _keep);
      intsItr++;
    }

  blacklisted.resize(noLocalParts);
  blacklisted.reset();
}

void ParticleManager::resizeAll()
{
  resizeAll(false);
}

///
/// resize a vectorial quantitiy
///
void ParticleManager::resize(matrixRefType _matrix, bool _keep)
{
  std::string matrixName = getName(_matrix);
  size_t newSize;

  if (isGhostVar(matrixName))
    {
      newSize = noLocalParts + noGhostParts;
    }
  else
    {
      newSize = noLocalParts;
    }

  if (_matrix.size1() != newSize)
    {
      _matrix.resize(newSize, _matrix.size2(), _keep);
    }
}

///
/// resize a scalar quantitiy
///
void ParticleManager::resize(valvectRefType _valvect, bool _keep)
{
  std::string vectName = getName(_valvect);
  size_t newSize;

  if (isGhostVar(vectName))
    {
      newSize = noLocalParts + noGhostParts;
    }
  else
    {
      newSize = noLocalParts;
    }

  if (_valvect.size() != newSize)
    {
      _valvect.resize(newSize, _keep);
    }
}

///
/// resize an integer quantitiy
///
void ParticleManager::resize(idvectRefType _idvect, bool _keep)
{
  std::string vectName = getName(_idvect);
  size_t newSize;

  if (isGhostVar(vectName))
    {
      newSize = noLocalParts + noGhostParts;
    }
  else
    {
      newSize = noLocalParts;
    }

  if (_idvect.size() != newSize)
    {
      _idvect.resize(newSize, _keep);
    }
}

///
/// return number of local or ghost particles
///
size_t ParticleManager::getNoLocalParts()
{
  return noLocalParts;
}

size_t ParticleManager::getNoGhostParts()
{
  return noGhostParts;
}

size_t ParticleManager::getNoTotalParts()
{
  return noLocalParts + noGhostParts;
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
  matrixPtrMap::iterator matrItr = usedVectors.find(_name);

  if (matrItr == usedVectors.end())
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
  valVectPtrMap::iterator vectItr = usedScalars.find(_name);

  if (vectItr == usedScalars.end())
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
  idVectPtrMap::iterator intItr = usedIntegers.find(_name);

  if (intItr == usedIntegers.end())
    {
      return NULL;
    }
  else
    {
      return intItr->second;
    }
}


///
/// register an auxilliary quantity, for example integration vars
///
void ParticleManager::regQuantity(matrixRefType _matrixRef,
                                  std::string _name)
{
  if ( _name.size() > 0 )
  {
    usedVectors[ _name ] = &_matrixRef;
  }
}

void ParticleManager::regQuantity(valvectRefType _valvectRef,
                                  std::string _name)
{
  if ( _name.size() > 0 )
  {
  usedScalars[ _name ] = &_valvectRef;
  }
}

void ParticleManager::regQuantity(idvectRefType _idvectRef,
                                  std::string _name)
{
  if ( _name.size() > 0 )
  {
  usedIntegers[ _name ] = &_idvectRef;
  }
}

///
/// unregister a quantity again, for example when an integrator is destructed
///
void ParticleManager::unRegQuantity(matrixRefType _matrRef)
{
  while ( getName(_matrRef).size() > 0 )
  {
    usedVectors.erase( getName(_matrRef) );
  }
}

void ParticleManager::unRegQuantity(valvectRefType _valvectRef)
{
  while ( getName(_valvectRef).size() > 0 )
  {
    usedScalars.erase( getName(_valvectRef) );
  }
}

void ParticleManager::unRegQuantity(idvectRefType _idvectRef)
{
  while ( getName(_idvectRef).size() > 0 )
  {
    usedIntegers.erase( getName(_idvectRef) );
  }
}

///
/// does the attribute exist?
///
bool ParticleManager::attrExists(std::string _name)
{
  return(attributes.find(_name) != attributes.end());
}
};

 #endif
