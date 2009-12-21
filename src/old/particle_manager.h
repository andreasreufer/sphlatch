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
/*enum VectIndices
{
   X, Y, Z
};*/
enum TensorIndices
{
   XX, XY, XZ, YY, YZ
};
}

class ParticleManager
{
public:
   typedef ParticleManager                         self_type;
   typedef ParticleManager&                        self_reference;
   typedef ParticleManager*                        self_pointer;

   typedef std::map<std::string, valvectPtrType>   valvectPtrMap;
   typedef std::map<std::string, matrixPtrType>    matrixPtrMap;
   typedef std::map<std::string, idvectPtrType>    idvectPtrMap;

   typedef std::set<std::string>                   stringSet;

   typedef std::map<std::string, std::string>      strStrMapType;

   static self_reference instance(void);

   ///
   /// methods to set up the vars needed
   ///
   void usePhaseSpace();
   void useBasicSPH();
   void useEnergy();
   void useTimedepEnergy();
   void useEntropy();
   void useTimedepEntropy();
   void useTimedepH();
   void useGravity();
   void useEccentricity();
   void useIntegratedRho();

   void useCost(void);

   void useAVMonaghan(void);
   void useAVBalsara(void);
   void useAVHernquist(void);
   void useAVTimedepAlpha(void);

   void useMaterials(void);
   void usePhase(void);
   void useTemperature(void);

   void useStress(void);
   void useDamage(void);

   void usePeakPress(void);

   void setNoParts(size_t _noParts, size_t _noGhostParts);
   void setNoParts(size_t _noParts);

   void addParts(size_t _addNoParts);

   void resizeAll(bool _keep);

   void resizeAll(void);

   void resize(matrixRefType _matrixRef, bool _keep);
   void resize(valvectRefType _valvectRef, bool _keep);
   void resize(idvectRefType _idvectRef, bool _keep);

   ///
   /// synonyms for var names
   ///
private:
   strStrMapType synonymes;

   ///
   /// add a synonym for a known name
   ///
public:
   void addSynonym(std::string _name, std::string _syn);
   std::string resolveSynonym(std::string _syn);

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
   void regQuantity(matrixRefType _matrixRef, std::string _name);
   void regQuantity(valvectRefType _valvectRef, std::string _name);
   void regQuantity(idvectRefType _idvectRef, std::string _name);

   ///
   /// unregister a quantity again, for example when an integrator
   /// is destructed
   ///
   void unRegQuantity(matrixRefType _matRef);
   void unRegQuantity(valvectRefType _valvectRef);
   void unRegQuantity(idvectRefType _idvectRef);

   ///
   /// vector quantities:
   ///
   matrixType pos,     /// position
              vel,     /// velocity
              acc,     /// acceleration
              rotv,    /// vorticity
              S, dSdt, /// traceless deviatoric stress tensor
              dR,      /// rotation rate tensor
              e;       /// strain rate tensor

   ///
   /// scalar quantities
   ///
   valvectType m,               /// mass
               h, dhdt,         /// smoothing length
               rho, drhodt,     /// density
               p,               /// pressurce
               u, dudt,         /// specific internal energy
               A, dAdt,         /// entropy
               alpha, dalphadt, /// artificial viscosity alpha
               divv,            /// velocity divergence
               mumax,           /// maximal Monaghan mu
               q,               /// AV parameter
               eps,             /// gravitational smoothing
               dtav,            /// AV timestep
               dt,              /// timestep
               cs,              /// speed of sound
               ecc,             /// eccentricity
               T,               /// temperature
               dam, ddamdt,     /// damage
               epsmin,          /// minimal strain
               acoef,           /// crack growth timescale
               mweib,           /// Weibull m exponent
               young,           /// Youngs modulus
               cost,            /// relative computational cost
               peakp;           /// peak pressure

   ///
   /// integers
   ///
   idvectType id,      /// ID
              noneigh, /// number of neighbours
              mat,     /// material
              phase,   /// phase
              noflaws; /// number of flaws

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
   quantsType usedQuants;

   matrixPtrMap  knownVects;
   valvectPtrMap knownScalars;
   idvectPtrMap  knownIntegers;

public:
   ///
   /// get/set used quantities
   ///
   quantsType getUsedQuants();
   void setUsedQuants(quantsType _quants);

   ///
   /// get known quants from a string list
   ///
   quantsType getKnownQuants(stringListType _strlist);

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
   stringSetType isGhostVarSet;
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
      return(*_instance);
   else
   {
      _instance = new ParticleManager;
      return(*_instance);
   }
}

ParticleManager::ParticleManager(void)
{
   using namespace boost::assign;
   ///
   /// known vectorial quantities with
   /// their respective second dimension
   /// size
   ///
   knownVects["pos"]  = &pos;
   knownVects["vel"]  = &vel;
   knownVects["acc"]  = &acc;
   knownVects["rotv"] = &rotv;
   knownVects["S"]    = &S;
   knownVects["dSdt"] = &dSdt;
   knownVects["dR"]   = &dR;
   knownVects["e"]    = &e;

   pos.resize(0, 3);
   vel.resize(0, 3);
   acc.resize(0, 3);
   rotv.resize(0, 3);
   S.resize(0, 5);
   dSdt.resize(0, 5);
   dR.resize(0, 5);
   e.resize(0, 5);

   ///
   /// vectorial quantities which are ghost
   ///
   isGhostVarSet += "pos", "vel", "rotv", "S";

   ///
   /// known integer quantities
   ///
   knownIntegers["id"]      = &id;
   knownIntegers["noneigh"] = &noneigh;
   knownIntegers["mat"]     = &mat;
   knownIntegers["phase"]   = &phase;
   knownIntegers["noflaws"] = &noflaws;

   ///
   /// integer   quantities which are ghost
   ///
   isGhostVarSet += "id";

   ///
   /// known scalar quantities
   ///
   knownScalars["m"]        = &m;
   knownScalars["h"]        = &h;
   knownScalars["dhdt"]     = &dhdt;
   knownScalars["rho"]      = &rho;
   knownScalars["drhodt"]   = &drhodt;
   knownScalars["p"]        = &p;
   knownScalars["u"]        = &u;
   knownScalars["dudt"]     = &dudt;
   knownScalars["A"]        = &A;
   knownScalars["dAdt"]     = &dAdt;
   knownScalars["alpha"]    = &alpha;
   knownScalars["dalphadt"] = &dalphadt;
   knownScalars["divv"]     = &divv;
   knownScalars["mumax"]    = &mumax;
   knownScalars["q"]        = &q;
   knownScalars["eps"]      = &eps;
   knownScalars["dtav"]     = &dtav;
   knownScalars["dt"]       = &dt;
   knownScalars["cs"]       = &cs;
   knownScalars["ecc"]      = &ecc;
   knownScalars["T"]        = &T;
   knownScalars["dam"]      = &dam;
   knownScalars["ddamdt"]   = &ddamdt;
   knownScalars["epsmin"]   = &epsmin;
   knownScalars["acoef"]    = &acoef;
   knownScalars["mweib"]    = &mweib;
   knownScalars["young"]    = &young;
   knownScalars["cost"]     = &cost;
   knownScalars["peakp"]    = &peakp;

   ///
   /// scalar quantities which are ghost
   ///
   isGhostVarSet += "m", "h", "rho", "p", "A",
   "alpha", "divv", "eps", "cs";

   ///
   /// set some initial values
   ///
   step = 0, substep = 0;

   noLocalParts = 0;
   noGhostParts = 0;
}

ParticleManager::~ParticleManager(void)
{ }

///
///
///
void ParticleManager::usePhaseSpace()
{
   using namespace boost::assign;
   usedQuants.vects += &pos, &vel;
}

///
/// everything for nbody
///
void ParticleManager::useGravity()
{
   usePhaseSpace();
   using namespace boost::assign;
   usedQuants.vects   += &acc;
   usedQuants.scalars += &m, &eps;
   usedQuants.ints    += &id;
}

void ParticleManager::useEccentricity()
{
   using namespace boost::assign;
   usedQuants.scalars += &ecc;
}

///
/// most basic SPH variables
///
void ParticleManager::useBasicSPH()
{
   usePhaseSpace();
   using namespace boost::assign;
   usedQuants.vects   += &acc;
   usedQuants.scalars += &m, &h, &rho, &p, &cs;
   usedQuants.ints    += &noneigh, &id;
}

///
/// use internal energy
///
void ParticleManager::useEnergy()
{
   using namespace boost::assign;
   usedQuants.scalars += &u;
}

///
/// time dependent internal energy
///
void ParticleManager::useTimedepEnergy()
{
   using namespace boost::assign;
   usedQuants.scalars += &dudt;
}

///
/// use entropy
///
void ParticleManager::useEntropy()
{
   using namespace boost::assign;
   usedQuants.scalars += &A;
}

///
/// time dependent entropy
///
void ParticleManager::useTimedepEntropy()
{
   using namespace boost::assign;
   usedQuants.scalars += &dAdt;
}

///
/// time dependent smoothing length
///
void ParticleManager::useTimedepH()
{
   using namespace boost::assign;
   usedQuants.scalars += &dhdt, &divv;
}

///
/// vars for Monaghan (1989) artificial viscosity
///
void ParticleManager::useAVMonaghan()
{
   using namespace boost::assign;
   usedQuants.scalars += &mumax;
}

///
/// time dependent alpha
///
void ParticleManager::useAVTimedepAlpha()
{
   using namespace boost::assign;
   usedQuants.scalars += &alpha, &dalphadt;
}

///
/// balsara switch
///
void ParticleManager::useAVBalsara()
{
   using namespace boost::assign;
   usedQuants.vects += &rotv;
}

///
/// vars for Hernquist&Katz artificial viscosity
///
void ParticleManager::useAVHernquist()
{
   using namespace boost::assign;
   usedQuants.scalars += &q;
}

///
/// the material variable
///
void ParticleManager::useMaterials()
{
   using namespace boost::assign;
   usedQuants.ints += &mat;
}

///
/// use the phase variable
///
void ParticleManager::usePhase()
{
   using namespace boost::assign;
   usedQuants.ints += &phase;
}

///
/// use the temperature variable
///
void ParticleManager::useTemperature()
{
   using namespace boost::assign;
   usedQuants.scalars += &T;
}

///
/// the density time derivative variable
///
void ParticleManager::useIntegratedRho()
{
   using namespace boost::assign;
   usedQuants.scalars += &drhodt;
}

///
/// the cost variable
///
void ParticleManager::useCost()
{
   using namespace boost::assign;
   usedQuants.scalars += &cost;
}

///
/// use variables for the stress tensor
///
void ParticleManager::useStress()
{
   using namespace boost::assign;
   usedQuants.vects += &S, &dSdt, &e, &dR;
}

///
/// use variables for the stress tensor
///
void ParticleManager::useDamage()
{
   using namespace boost::assign;
   usedQuants.scalars += &dam, &ddamdt, &epsmin, &acoef, &mweib, &young;
   usedQuants.ints    += &noflaws;
}

///
/// use variables for the stress tensor
///
void ParticleManager::usePeakPress()
{
   using namespace boost::assign;
   usedQuants.scalars += &peakp;
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
   matrixPtrSetType::iterator vectorsItr = usedQuants.vects.begin();

   while (vectorsItr != usedQuants.vects.end())
   {
      resize(**vectorsItr, _keep);
      vectorsItr++;
   }

   valvectPtrSetType::iterator scalarsItr = usedQuants.scalars.begin();
   while (scalarsItr != usedQuants.scalars.end())
   {
      resize(**scalarsItr, _keep);
      scalarsItr++;
   }

   idvectPtrSetType::iterator intsItr = usedQuants.ints.begin();
   while (intsItr != usedQuants.ints.end())
   {
      resize(**intsItr, _keep);
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
   size_t      newSize;

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
   size_t      newSize;

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
   size_t      newSize;

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
   return(noLocalParts);
}

size_t ParticleManager::getNoGhostParts()
{
   return(noGhostParts);
}

size_t ParticleManager::getNoTotalParts()
{
   return(noLocalParts + noGhostParts);
}

///
/// get/set used quantitites
///
quantsType ParticleManager::getUsedQuants()
{
   return(usedQuants);
}

void ParticleManager::setUsedQuants(quantsType _quants)
{
   usedQuants = _quants;
}

///
/// get known quants from a string list
///
quantsType ParticleManager::getKnownQuants(stringListType _strlist)
{
   quantsType retQuants;

   stringListType::iterator strItr = _strlist.begin();

   while (strItr != _strlist.end())
   {
      matrixPtrMap::iterator vectItr = knownVects.find(*strItr);
      if (vectItr != knownVects.end())
         retQuants.vects.insert(vectItr->second);

      valvectPtrMap::iterator scalItr = knownScalars.find(*strItr);
      if (scalItr != knownScalars.end())
         retQuants.scalars.insert(scalItr->second);

      idvectPtrMap::iterator intItr = knownIntegers.find(*strItr);
      if (intItr != knownIntegers.end())
         retQuants.ints.insert(intItr->second);

      strItr++;
   }

   return(retQuants);
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
///
///
void ParticleManager::addSynonym(std::string _name, std::string _syn)
{
   bool isKnown = false;
   valvectPtrMap::iterator vectItr = knownScalars.find(_name);

   if (vectItr != knownScalars.end())
      isKnown = true;

   matrixPtrMap::iterator matrItr = knownVects.find(_name);
   if (matrItr != knownVects.end())
      isKnown = true;

   idvectPtrMap::iterator idItr = knownIntegers.find(_name);
   if (idItr != knownIntegers.end())
      isKnown = true;

   if (isKnown == true)
      synonymes[_syn] = _name;
}

std::string ParticleManager::resolveSynonym(std::string _syn)
{
   strStrMapType::iterator synItr = synonymes.find(_syn);

   if (synItr != synonymes.end())
      return(synItr->second);
   else
      return(_syn);
}

///
/// get a vectorial quantitiy references name
///
std::string ParticleManager::getName(matrixRefType _matrixRef)
{
   matrixPtrSetType::iterator usedItr = usedQuants.vects.find(&_matrixRef);

   if (usedItr != usedQuants.vects.end())
   {
      matrixPtrMap::iterator vectsItr = knownVects.begin();
      while (vectsItr != knownVects.end())
      {
         if (vectsItr->second == &_matrixRef)
            return(vectsItr->first);

         vectsItr++;
      }
      return("");
   }
   else
      return("");
}

///
/// get a scalar quantitiy references name
///
std::string ParticleManager::getName(valvectRefType _valvectRef)
{
   valvectPtrSetType::iterator usedItr =
      usedQuants.scalars.find(&_valvectRef);

   if (usedItr != usedQuants.scalars.end())
   {
      valvectPtrMap::iterator scalItr = knownScalars.begin();
      while (scalItr != knownScalars.end())
      {
         if (scalItr->second == &_valvectRef)
            return(scalItr->first);

         scalItr++;
      }
      return("");
   }
   else
      return("");
}

///
/// get an integer quantitiy references name
///
std::string ParticleManager::getName(idvectRefType _idvectRef)
{
   idvectPtrSetType::iterator usedItr = usedQuants.ints.find(&_idvectRef);

   if (usedItr != usedQuants.ints.end())
   {
      idvectPtrMap::iterator intItr = knownIntegers.begin();
      while (intItr != knownIntegers.end())
      {
         if (intItr->second == &_idvectRef)
            return(intItr->first);

         intItr++;
      }
      return("");
   }
   else
      return("");
}

///
/// get a vectorial quantity pointer by name
///
matrixPtrType ParticleManager::getVectRef(std::string _name)
{
   _name = resolveSynonym(_name);
   matrixPtrMap::iterator matrItr = knownVects.find(_name);

   if (matrItr == knownVects.end())
      return(NULL);
   else
   {
      if (usedQuants.vects.count(matrItr->second) > 0)
         return(matrItr->second);
      else
         return(NULL);
   }
}

///
/// get a scalar quantity pointer by name
///
valvectPtrType ParticleManager::getScalarRef(std::string _name)
{
   _name = resolveSynonym(_name);
   valvectPtrMap::iterator vectItr = knownScalars.find(_name);

   if (vectItr == knownScalars.end())
      return(NULL);
   else
   {
      if (usedQuants.scalars.count(vectItr->second) > 0)
         return(vectItr->second);
      else
         return(NULL);
   }
}

///
/// get an integer quantity pointer by name
///
idvectPtrType ParticleManager::getIdRef(std::string _name)
{
   _name = resolveSynonym(_name);
   idvectPtrMap::iterator intItr = knownIntegers.find(_name);

   if (intItr == knownIntegers.end())
      return(NULL);
   else
   {
      if (usedQuants.ints.count(intItr->second) > 0)
         return(intItr->second);
      else
         return(NULL);
   }
}

///
/// register an auxilliary quantity, for example integration vars
///
void ParticleManager::regQuantity(matrixRefType _matrixRef,
                                  std::string   _name)
{
   knownVects[_name] = &_matrixRef;
   usedQuants.vects.insert(&_matrixRef);
}

void ParticleManager::regQuantity(valvectRefType _valvectRef,
                                  std::string    _name)
{
   knownScalars[_name] = &_valvectRef;
   usedQuants.scalars.insert(&_valvectRef);
}

void ParticleManager::regQuantity(idvectRefType _idvectRef,
                                  std::string   _name)
{
   knownIntegers[_name] = &_idvectRef;
   usedQuants.ints.insert(&_idvectRef);
}

///
/// unregister a quantity again, for example when an integrator is destructed
///
void ParticleManager::unRegQuantity(matrixRefType _matrRef)
{
   knownVects.erase( getName( _matrRef ) );
   usedQuants.vects.erase(&_matrRef);
}

void ParticleManager::unRegQuantity(valvectRefType _valvectRef)
{
   knownScalars.erase( getName( _valvectRef ) );
   usedQuants.scalars.erase(&_valvectRef);
}

void ParticleManager::unRegQuantity(idvectRefType _idvectRef)
{
   knownIntegers.erase( getName( _idvectRef ) );
   usedQuants.ints.erase(&_idvectRef);
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
