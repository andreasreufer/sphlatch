#ifndef SPHLATCH_COSTZONE_H
#define SPHLATCH_COSTZONE_H

#include <iostream>
#include <vector>

#include "communicationmanager.h"
#include "memorymanager.h"

namespace sphlatch
{

class CostZone
{
public:

typedef CostZone self_type;
typedef CostZone& self_reference;
typedef CostZone* self_pointer;

typedef sphlatch::CommunicationManager commManagerType;
typedef sphlatch::MemoryManager memoryManagerType;

private:

static self_pointer _instance;
commManagerType& CommManager;
memoryManagerType& MemManager;

public:

static self_reference instance(void);

domainPartsIndexRefType createDomainIndexVector(void);
domainPartsIndexRefType createDomainGhostIndexVector(void);

domainPartsIndexType domainGhostIndexVector;

valueType getSidelength();
valvectType getCenter();
size_t getDepth();

void resize(const size_t size);

protected:

CostZone();
~CostZone(void);


private:
valueType xmin, ymin, zmin, xmax, ymax, zmax, sidelength;

};

CostZone::self_pointer CostZone::_instance = NULL;

CostZone::self_reference CostZone::instance(void)
{
  if (_instance != NULL)
    return *_instance;
  else
    {
      _instance = new CostZone;
      return *_instance;
    }
}

CostZone::CostZone(void) :
  CommManager(commManagerType::instance()),
  MemManager(memoryManagerType::instance())
{
}

CostZone::~CostZone(void)
{
}

domainPartsIndexRefType CostZone::createDomainIndexVector(void)
{
  return domainGhostIndexVector;
};

domainPartsIndexRefType CostZone::createDomainGhostIndexVector(void)
{
  return domainGhostIndexVector;
}

valvectType CostZone::getCenter(void)
{
  valvectType retvect(3);
  retvect(0) = ( xmin + xmax ) / 2;
  retvect(1) = ( ymin + ymax ) / 2;
  retvect(2) = ( zmin + zmax ) / 2;
  return retvect;
}

valueType CostZone::getSidelength(void)
{
  return sidelength;
}


};

#endif
