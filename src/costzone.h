#ifndef SPHLATCH_COSTZONE_H
#define SPHLATCH_COSTZONE_H

#include <vector>
#include <limits>

#include "communicationmanager.h"
#include "memorymanager.h"
#include "spacefillingcurve.h"

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

//#define SPHLATCH_CARTESIAN_XYZ
//#define SPHLATCH_CARTESIAN_YZX
//#define SPHLATCH_CARTESIAN_ZXY
#define SPHLATCH_HILBERT3D

#ifdef SPHLATCH_CARTESIAN_XYZ
typedef SpaceFillingCurve<CartesianXYZ> sfcurveType;
#endif
#ifdef SPHLATCH_CARTESIAN_YZX
typedef SpaceFillingCurve<CartesianYZX> sfcurveType;
#endif
#ifdef SPHLATCH_CARTESIAN_ZXY
typedef SpaceFillingCurve<CartesianZXY> sfcurveType;
#endif
#ifdef SPHLATCH_HILBERT2D_XY

#endif
#ifdef SPHLATCH_HILBERT2D_XZ

#endif
#ifdef SPHLATCH_HILBERT2D_YZ

#endif
#ifdef SPHLATCH_HILBERT3D
typedef SpaceFillingCurve<Hilbert3D> sfcurveType;
#endif

private:

void centerOfTheUniverse();

static self_pointer _instance;
commManagerType& CommManager;
memoryManagerType& MemManager;
sfcurveType spaceCurve;

public:

static self_reference instance(void);

domainPartsIndexRefType createDomainPartsIndex(void);
domainPartsIndexRefType createDomainGhostIndex(void);

static domainPartsIndexType domainPartsIndex;
static domainPartsIndexType domainGhostIndex;

valueType getSidelength();
valvectType getCenter();
size_t getDepth();

void resize(const size_t size);

protected:

CostZone();
~CostZone(void);


private:
static domainPartsIndexType costzoneCells;
static countsVectType partsPerCell;
valueType xCenter, yCenter, zCenter, sidelength;
size_t depth, noCells1D, noCells2D, noCells3D;
size_t myFirstWalkIndex, myLastWalkIndex;
};

CostZone::self_pointer CostZone::_instance = NULL;

domainPartsIndexType CostZone::domainPartsIndex;
domainPartsIndexType CostZone::domainGhostIndex;
domainPartsIndexType CostZone::costzoneCells;
countsVectType CostZone::partsPerCell;

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
///
/// standard costzone depth
///
  //resize(3);
  resize(5);
  //resize(6);
}

CostZone::~CostZone(void)
{
}

void CostZone::resize(size_t _depth)
{
  depth = _depth;
  noCells1D = (1 << depth); /// 2^depth
  noCells2D = noCells1D * noCells1D;
  noCells3D = noCells2D * noCells1D;

  ///
  /// resize various containers
  ///
  costzoneCells.resize(noCells3D);
  partsPerCell.resize(noCells3D);

  spaceCurve.init(depth);
}

///
/// returns a vector of index lists (vector) of the particles
/// belonging to each domain
/// note that the domain numbers do not have to coincidence
/// with the MPI rank
///
domainPartsIndexRefType CostZone::createDomainPartsIndex(void)
{
  matrixRefType Data(MemManager.Data);

  const size_t noDomains = CommManager.getNoDomains();
  const size_t myDomain = CommManager.getMyDomain();
  const size_t noParts = Data.size1();

  domainPartsIndex.resize(noDomains);

  centerOfTheUniverse();

  ///
  /// prepare the costzone cells
  ///
  for (size_t i = 0; i < noCells3D; i++)
    {
      costzoneCells[i].resize(0);
      partsPerCell[i] = 0;
    }

  ///
  /// define the position corresponding to index 0
  ///
  const valueType xMin = xCenter - 0.5 * sidelength;
  const valueType yMin = yCenter - 0.5 * sidelength;
  const valueType zMin = zCenter - 0.5 * sidelength;
  ///
  /// this factor translates distance to cell index
  ///
  const valueType lengthToIndex = static_cast<valueType>(noCells1D - 1) / sidelength;

  ///
  /// put the local particles in the corresponding costzone cells
  ///
  for (size_t i = 0; i < noParts; i++)
    {
      const size_t xIndex = lrint((Data(i, X) - xMin) * lengthToIndex);
      const size_t yIndex = lrint((Data(i, Y) - yMin) * lengthToIndex);
      const size_t zIndex = lrint((Data(i, Z) - zMin) * lengthToIndex);

      const size_t cartIndex = xIndex + yIndex * noCells1D + zIndex * noCells2D;

      costzoneCells[cartIndex].push_back(i);
      partsPerCell[cartIndex] += 1;
    }

  ///
  /// globally sum up partsPerCell and get
  /// the global number of particles
  ///
  countsVectType locPartPerCell = partsPerCell;

  CommManager.sumUpCounts(partsPerCell);
  size_t noGlobParts = 0;
  for (size_t i = 0; i < noCells3D; i++)
    {
      noGlobParts += partsPerCell[i];
    }

  ///
  /// determine the amount of particles each domain should have
  ///
  countsVectType domainSize;
  domainSize.resize(noDomains);
  for (size_t i = 0; i < noDomains; i++)
    {
      domainSize[i] = lrint(static_cast<double>(noGlobParts) /
                            static_cast<double>(noDomains));
    }

  ///
  /// now let's walk along the spacefilling curve
  ///
  size_t curDomain = 0;
  size_t noDistrParts = 0;
  size_t noCurDomainParts = 0;
  for (size_t i = 0; i < noCells3D; i++)
    {
      const size_t cartIndex = spaceCurve.curveIndexToCartIndex(i);

      ///
      /// this algorithm is not ideal, revise ...
      /// if ( addition of current cell completes the distribution OD
      ///      the current domain is saturated )
      ///    AND the current domain is not the last one
      ///    AND the current domain already has particles
      /// then start filling the next domain
      ///
      if ((noDistrParts + partsPerCell[cartIndex] == noGlobParts ||
           noCurDomainParts > static_cast<size_t>(domainSize[curDomain]))
          && curDomain < (noDomains - 1)
          && noCurDomainParts > 0)
        {
          curDomain++;
          noCurDomainParts = 0;
        }

      const size_t noCurParts = partsPerCell[cartIndex];
      noCurDomainParts += noCurParts;
      noDistrParts += noCurParts;

      ///
      /// use vector iterators
      for (size_t i = 0; i < costzoneCells[cartIndex].size(); i++)
        {
          domainPartsIndex[curDomain].push_back(costzoneCells[cartIndex][i]);
        }
    }

  return domainPartsIndex;
};

///
/// returns a vector of index lists (vector) of the particles
/// belonging to each domain
/// note that the domain numbers do not have to coincidence
/// with the MPI rank
///
domainPartsIndexRefType CostZone::createDomainGhostIndex(void)
{
  return domainGhostIndex;;
}

///
/// calculates the center of the universe
/// yes, it can be done :-)
///
void CostZone::centerOfTheUniverse(void)
{
  matrixRefType Data(MemManager.Data);
  const size_t noParts = Data.size1();

  valueType xMin = std::numeric_limits<valueType>::max();
  valueType yMin = std::numeric_limits<valueType>::max();
  valueType zMin = std::numeric_limits<valueType>::max();

  valueType xMax = std::numeric_limits<valueType>::min();
  valueType yMax = std::numeric_limits<valueType>::min();
  valueType zMax = std::numeric_limits<valueType>::min();

  for (size_t i = 0; i < noParts; i++)
    {
      xMin = Data(i, X) < xMin ? Data(i, X) : xMin;
      xMax = Data(i, X) > xMax ? Data(i, X) : xMax;

      yMin = Data(i, Y) < yMin ? Data(i, Y) : yMin;
      yMax = Data(i, Y) > yMax ? Data(i, Y) : yMax;

      zMin = Data(i, Z) < zMin ? Data(i, Z) : zMin;
      zMax = Data(i, Z) > zMax ? Data(i, Z) : zMax;
    }

  CommManager.min(xMin);
  CommManager.max(xMax);

  CommManager.min(yMin);
  CommManager.max(yMax);

  CommManager.min(zMin);
  CommManager.max(zMax);

  sidelength = std::max(std::max(xMax - xMin, yMax - yMin), zMax - zMin);

  xCenter = (xMin + xMax) / 2.;
  yCenter = (yMin + yMax) / 2.;
  zCenter = (zMin + zMax) / 2.;
}

valvectType CostZone::getCenter(void)
{
  valvectType retvect(3);

  retvect(0) = xCenter;
  retvect(1) = yCenter;
  retvect(2) = zCenter;
  return retvect;
}

valueType CostZone::getSidelength(void)
{
  return sidelength;
}
};

#endif
