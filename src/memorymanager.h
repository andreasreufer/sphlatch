#ifndef SPHLATCH_MEMORYMANAGER_H
#define SPHLATCH_MEMORYMANAGER_H

#include <map>
#include <set>
#include <boost/assign/std/vector.hpp>

#include "typedefs.h"
#include "particle.h"

namespace sphlatch
{
class MemoryManager
{
public:

typedef MemoryManager self_type;
typedef MemoryManager& self_reference;
typedef MemoryManager* self_pointer;

typedef std::map<std::string, valueType> strValMapType;

matrixType Data, GData;

std::string simName;

/** PopParticle pops the particle with
*  a given row index or particle ID */
valvectType popParticle(valueType _particleID);
valvectType popParticle(size_t _particleIndex);

/** PopParticles pops particles with
 *  given row indeces or particle ID */
matrixType popParticles(valueType _particleID);
matrixType popParticles(idvectType _particleIndices);

/** DelParticle(s) deletes one or more
 *  particle with given row indices */
void delParticle(size_t _particleIndex);
void delParticles(std::vector<size_t> _particleIndices);

/** PushParticle() pushes a particle
 *  to the end of Data and returns the row index */
size_t pushParticle(valvectType _particle);

bool saveParameter(std::string _paramName, valueType _paramValue, bool _queueParam);
valueType loadParameter(std::string _paramName);
std::set<std::string> dumpParameter(void);

static self_reference instance(void);

protected:

MemoryManager(void);
~MemoryManager(void);

private:

static self_pointer _instance;
strValMapType paramMap;
std::set<std::string> paramSave;
};


MemoryManager::self_pointer MemoryManager::_instance = NULL;

MemoryManager::self_reference MemoryManager::instance(void)
{
  if (_instance != NULL)
    return *_instance;
  else
    {
      _instance = new MemoryManager;
      return *_instance;
    }
}

MemoryManager::MemoryManager(void)
{
}

MemoryManager::~MemoryManager(void)
{
}

valueType MemoryManager::loadParameter(std::string _paramName)
{
  if (paramMap.count(_paramName) > 0)
    {
      return paramMap[_paramName];
    }
  else
    {
      return std::numeric_limits<valueType>::quiet_NaN();
    }
}

bool MemoryManager::saveParameter(std::string _paramName,
                                  valueType _paramValue,
                                  bool _queueParam)
{
  paramMap[_paramName] = _paramValue;
  if (_queueParam)
    {
      paramSave.insert(_paramName);
    }
  return true;
}

std::set<std::string> MemoryManager::dumpParameter()
{
  return paramSave;
}

valvectType MemoryManager::popParticle(valueType _particleID)
{
  const size_t noParts = Data.size1();
  bool PartFound = false;
  size_t i = 0;
  valvectType PopParticle;

  if (noParts > 0)
    {
      while (!PartFound && i < noParts)
        {
          if (Data(i, ID) == _particleID)                   // not clean
            {
              PartFound = true;

              particleRowType MatchedParticle(Data, i);
              particleRowType LastParticle(Data, noParts - 1);
              PopParticle = MatchedParticle;

              LastParticle.swap(MatchedParticle);
              Data.resize(noParts - 1, Data.size2(), true);
            }
          i++;
        }
    }

  return PopParticle;
}

valvectType MemoryManager::popParticle(size_t _particleIndex)
{
  const size_t noParts = Data.size1();
  valvectType PopParticle;

  if (_particleIndex < noParts)
    {
      particleRowType MatchedParticle(Data, _particleIndex);
      particleRowType LastParticle(Data, noParts - 1);
      PopParticle = MatchedParticle;

      LastParticle.swap(MatchedParticle);
      Data.resize(noParts - 1, Data.size2(), true);
    }

  return PopParticle;
}

matrixType MemoryManager::popParticles(idvectType _particleIndices)
{
  const size_t noParts = Data.size1();
  const size_t NoVars = Data.size2();
  const size_t NoPopParticles = _particleIndices.size();

  matrixType PoppedData;

  PoppedData.resize(NoPopParticles, NoVars);

  size_t ParticlesPopped = 0;

  std::sort(_particleIndices.begin(), _particleIndices.end(),
            std::greater<size_t>());

  for (size_t i = 0; i < NoPopParticles; i++)
    {
      const size_t ParticleIndex = _particleIndices[i];
      if (ParticleIndex < noParts)
        {
          particleRowType MatchedParticle(Data, ParticleIndex);
          particleRowType(PoppedData, ParticlesPopped) = MatchedParticle;
          particleRowType LastParticle(Data, noParts - 1 - ParticlesPopped);
          LastParticle.swap(MatchedParticle);
          ParticlesPopped++;
        }
    }

  Data.resize(noParts - ParticlesPopped, NoVars, true);
  PoppedData.resize(ParticlesPopped, NoVars, true);

  return PoppedData;
}


matrixType MemoryManager::popParticles(valueType _particleID)
{
  using namespace boost::assign;

  const size_t noParts = Data.size1();
  const size_t NoVars = Data.size2();

  std::vector<size_t> ParticleIndices;

  for (int i = noParts - 1; i >= 0; --i)
    {
      if (lrint(Data(i, ID)) == lrint(_particleID))
        {
          ParticleIndices += i;
        }
    }

  const size_t NoPopParticles = ParticleIndices.size();

  matrixType PoppedData;
  PoppedData.resize(NoPopParticles, NoVars);

  for (size_t i = 0; i < NoPopParticles; i++)
    {
      const size_t ParticleIndex = ParticleIndices[i];

      particleRowType MatchedParticle(Data, ParticleIndex);
      particleRowType(PoppedData, i) = MatchedParticle;
      particleRowType LastParticle(Data, noParts - 1 - i);
      LastParticle.swap(MatchedParticle);
    }

  Data.resize(noParts - NoPopParticles, NoVars, true);
  PoppedData.resize(NoPopParticles, NoVars, true);

  return PoppedData;
}


void MemoryManager::delParticle(size_t _particleIndex)
{
  const size_t noParts = Data.size1();

  if (_particleIndex < noParts)
    {
      particleRowType MatchedParticle(Data, _particleIndex);
      particleRowType LastParticle(Data, noParts - 1);
      LastParticle.swap(MatchedParticle);
      Data.resize(noParts - 1, Data.size2(), true);
    }
}

void MemoryManager::delParticles(std::vector<size_t> _particleIndices)
{
  const size_t noParts = Data.size1();
  const size_t NoDelParticles = _particleIndices.size();
  size_t ParticlesDeleted = 0;

  std::sort(_particleIndices.begin(),
            _particleIndices.end(),
            std::greater<size_t>());

  for (size_t i = 0; i < NoDelParticles; i++)
    {
      const size_t ParticleIndex = _particleIndices[i];
      if (ParticleIndex < noParts)
        {
          particleRowType MatchedParticle(Data, ParticleIndex);
          particleRowType LastParticle(Data, noParts - 1 - ParticlesDeleted);
          LastParticle.swap(MatchedParticle);
          ParticlesDeleted++;
        }
    }
  Data.resize(noParts - ParticlesDeleted, Data.size2(), true);
}

size_t MemoryManager::pushParticle(valvectType _particle)
{
  const size_t InsertionIndex = Data.size1();
  const size_t PARTSIZE = Data.size2();

  Data.resize(InsertionIndex + 1, PARTSIZE, true);
  _particle.resize(PARTSIZE);

  particleRowType(Data, InsertionIndex) = _particle;

  return InsertionIndex;
}
};

#endif
