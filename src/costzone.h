#ifndef SPHLATCH_COSTZONE_H
#define SPHLATCH_COSTZONE_H

#include <vector>
#include <set>
#include <limits>
#include <cassert>

#include "typedefs.h"
#ifdef SPHLATCH_PARALLEL
#include "communication_manager.h"
#endif
#include "particle_manager.h"
#include "spacefillingcurve.h"

using namespace sphlatch::vectindices;

namespace sphlatch
{
class CostZone
{
public:

typedef CostZone self_type;
typedef CostZone& self_reference;
typedef CostZone* self_pointer;

#ifdef SPHLATCH_PARALLEL
typedef sphlatch::CommunicationManager commManagerType;
#endif
typedef sphlatch::ParticleManager partManagerType;

//#define SPHLATCH_CARTESIAN_XYZ
//#define SPHLATCH_CARTESIAN_YZX
//#define SPHLATCH_CARTESIAN_ZXY
//#define SPHLATCH_HILBERT3D

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
void fillCostzoneCells();

static self_pointer _instance;
#ifdef SPHLATCH_PARALLEL
commManagerType& CommManager;
#endif
partManagerType& PartManager;
sfcurveType spaceCurve;

public:

static self_reference instance(void);

domainPartsIndexRefType createDomainPartsIndex(void);
domainPartsIndexRefType createDomainGhostIndex(void);

static domainPartsIndexType domainPartsIndex;
static domainPartsIndexType domainGhostIndex;

valueType getSidelength();
valueType getAtomicLength();
valvectType getCenter();
size_t getDepth();
size_t getNoGhosts();

void resize(const size_t size);

protected:

CostZone();
~CostZone(void);


private:
static domainPartsIndexType costzoneCells;
static countsVectType partsPerCell;
static countsVectType domainMap;
valueType xCenter, yCenter, zCenter, sidelength;
size_t depth, noCells1D, noCells2D, noCells3D, noGhosts;
size_t myFirstWalkIndex, myLastWalkIndex;
};

CostZone::self_pointer CostZone::_instance = NULL;

domainPartsIndexType CostZone::domainPartsIndex;
domainPartsIndexType CostZone::domainGhostIndex;
domainPartsIndexType CostZone::costzoneCells;
countsVectType CostZone::partsPerCell;
countsVectType CostZone::domainMap;

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
#ifdef SPHLATCH_PARALLEL
  CommManager(commManagerType::instance()),
#endif
  PartManager(partManagerType::instance())
{
  ///
  /// standard costzone depth
  ///
  //const size_t defaultDepth = 1;
  //const size_t defaultDepth = 2;
  //const size_t defaultDepth = 3;
  //const size_t defaultDepth = 4;

  const size_t defaultDepth = 5;

  resize(defaultDepth);

#ifdef SPHLATCH_PARALLEL
  const size_t noDomains = CommManager.getNoDomains();
#else
  const size_t noDomains = 1;
#endif
  domainPartsIndex.resize(noDomains);
  domainGhostIndex.resize(noDomains);

  noGhosts = 0;
}

CostZone::~CostZone(void)
{
}

void CostZone::resize(size_t _depth)
{
  depth = _depth;
  noCells1D = (1 << depth);   /// 2^depth
  noCells2D = noCells1D * noCells1D;
  noCells3D = noCells2D * noCells1D;

  ///
  /// resize various containers
  ///
  costzoneCells.resize(noCells3D);
  partsPerCell.resize(noCells3D);
  domainMap.resize(noCells3D);

  spaceCurve.init(depth);
}

///
/// returns a vector of domains, containing all the indices of the
/// locally residing particles, which are assigned to the corresponding
/// domain
///
/// note that the domain numbers do not have to coincidence
/// with the MPI rank
///
domainPartsIndexRefType CostZone::createDomainPartsIndex(void)
{
#ifdef SPHLATCH_PARALLEL
  const size_t noDomains = CommManager.getNoDomains();
  const size_t myDomain = CommManager.getMyDomain();
#else
  const size_t noDomains = 1;
  const size_t myDomain = 0;
#endif

  for (size_t i = 0; i < noDomains; i++)
    {
      domainPartsIndex[i].resize(0);
    }

  ///
  /// determine the size of the universe
  /// and fill the costzone cells
  ///
  centerOfTheUniverse();
  fillCostzoneCells();

  ///
  /// globally sum up partsPerCell and get
  /// the global number of particles
  ///
  countsVectType globPartsPerCell = partsPerCell;

#ifdef SPHLATCH_PARALLEL
  CommManager.sumUpCounts(globPartsPerCell);
#endif
  size_t noGlobParts = 0;
  for (size_t i = 0; i < noCells3D; i++)
    {
      noGlobParts += globPartsPerCell[i];
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

  myFirstWalkIndex = noCells3D;
  myLastWalkIndex = noCells3D - 1;

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
      if ((noDistrParts + globPartsPerCell[cartIndex] == noGlobParts ||
           noCurDomainParts > static_cast<size_t>(domainSize[curDomain]))
          && curDomain < (noDomains - 1)
          && noCurDomainParts > 0)
        {
          if (myDomain == curDomain)
            {
              myLastWalkIndex = i - 1;
            }
          curDomain++;
          noCurDomainParts = 0;
          if (myDomain == curDomain)
            {
              myFirstWalkIndex = i;
            }
        }

      const size_t noCurParts = globPartsPerCell[cartIndex];
      noCurDomainParts += noCurParts;
      noDistrParts += noCurParts;

      domainMap[cartIndex] = curDomain;

      partsIndexVectType::const_iterator partsItr = (costzoneCells[cartIndex]).begin();
      const partsIndexVectType::const_iterator lastPart = (costzoneCells[cartIndex]).end();
      while (partsItr != lastPart)
        {
          domainPartsIndex[curDomain].push_back(*partsItr);
          partsItr++;
        }
    }

  ///
  /// set myFirstWalkIndex for domain 0 and
  ///
  if (myDomain == 0)
    {
      myFirstWalkIndex = 0;
    }

  ///
  /// determine number of expected ghosts
  ///
  countsVectType noGlobGhosts(noDomains);

  for (size_t i = myFirstWalkIndex; i <= myLastWalkIndex; i++)
    {
      const size_t cartIndex = spaceCurve.curveIndexToCartIndex(i);

      ///
      /// only check for neighbouring non-local cells, if the
      /// current cell actually has particles in it
      ///
      if (globPartsPerCell[cartIndex] > 0)
        {
          const size_t xCen = cartIndex % noCells1D;
          const size_t yCen = ((cartIndex - xCen) % noCells2D) / noCells1D;
          const size_t zCen = (cartIndex - xCen - yCen * noCells1D) / noCells2D;

          const size_t xMin = std::max(static_cast<int>(xCen - 1), 0);
          const size_t yMin = std::max(static_cast<int>(yCen - 1), 0);
          const size_t zMin = std::max(static_cast<int>(zCen - 1), 0);

          const size_t xMax = std::min(xCen + 1, noCells1D - 1);
          const size_t yMax = std::min(yCen + 1, noCells1D - 1);
          const size_t zMax = std::min(zCen + 1, noCells1D - 1);

          /// make this set static?
          std::set<size_t> domainSet;
          for (size_t xIndex = xMin; xIndex <= xMax; xIndex++)
            {
              for (size_t yIndex = yMin; yIndex <= yMax; yIndex++)
                {
                  for (size_t zIndex = zMin; zIndex <= zMax; zIndex++)
                    {
                      const size_t remCartIndex = xIndex + yIndex * noCells1D +
                                                  zIndex * noCells2D;
                      const size_t remDomain = domainMap[remCartIndex];
                      if (remDomain != myDomain)
                        {
                          domainSet.insert(remDomain);
                        }
                    }
                }
            }

          std::set<size_t>::const_iterator domainItr = domainSet.begin();
          const std::set<size_t>::const_iterator lastDomain = domainSet.end();

          while (domainItr != lastDomain)
            {
              noGlobGhosts[*domainItr] += globPartsPerCell[cartIndex];
              domainItr++;
            }
        }
    }
#ifdef SPHLATCH_PARALLEL
  CommManager.sumUpCounts(noGlobGhosts);
#endif
  noGhosts = noGlobGhosts[myDomain];

  return domainPartsIndex;
}
;

///
/// returns a vector of local particles indices, which
/// are ghosts to other domains
/// this function assumes, that createDomainPartsIndex() has been
/// executed beforehand and therefore only particles assigned to
/// the local domain actually reside on the local domain. this also
/// guarantees, that the domainMap is correctly filled
///
/// note that the domain numbers do not have to coincidence
/// with the MPI rank
///
domainPartsIndexRefType CostZone::createDomainGhostIndex(void)
{
#ifdef SPHLATCH_PARALLEL
  const size_t noDomains = CommManager.getNoDomains();
  const size_t myDomain = CommManager.getMyDomain();
#else
  const size_t noDomains = 1;
  const size_t myDomain = 0;
#endif

  for (size_t i = 0; i < noDomains; i++)
    {
      domainGhostIndex[i].resize(0);
    }

  ///
  /// fill again the costzone cells
  ///
  fillCostzoneCells();

  for (size_t i = myFirstWalkIndex; i <= myLastWalkIndex; i++)
    {
      const size_t cartIndex = spaceCurve.curveIndexToCartIndex(i);

      ///
      /// only check for neighbouring non-local cells, if the
      /// current cell actually has particles in it
      ///
      if (partsPerCell[cartIndex] > 0)
        {
          const size_t xCen = cartIndex % noCells1D;
          const size_t yCen = ((cartIndex - xCen) % noCells2D) / noCells1D;
          const size_t zCen = (cartIndex - xCen - yCen * noCells1D) / noCells2D;

          const size_t xMin = std::max(static_cast<int>(xCen - 1), 0);
          const size_t yMin = std::max(static_cast<int>(yCen - 1), 0);
          const size_t zMin = std::max(static_cast<int>(zCen - 1), 0);

          const size_t xMax = std::min(xCen + 1, noCells1D - 1);
          const size_t yMax = std::min(yCen + 1, noCells1D - 1);
          const size_t zMax = std::min(zCen + 1, noCells1D - 1);

          /// make this set static?
          std::set<size_t> domainSet;
          for (size_t xIndex = xMin; xIndex <= xMax; xIndex++)
            {
              for (size_t yIndex = yMin; yIndex <= yMax; yIndex++)
                {
                  for (size_t zIndex = zMin; zIndex <= zMax; zIndex++)
                    {
                      const size_t remCartIndex = xIndex + yIndex * noCells1D +
                                                  zIndex * noCells2D;
                      const size_t remDomain = domainMap[remCartIndex];
                      if (remDomain != myDomain)
                        {
                          domainSet.insert(remDomain);
                        }
                    }
                }
            }

          partsIndexVectType::const_iterator partsItr;
          const partsIndexVectType::const_iterator firstPart =
            (costzoneCells[cartIndex]).begin();
          const partsIndexVectType::const_iterator lastPart =
            (costzoneCells[cartIndex]).end();

          std::set<size_t>::const_iterator domainItr = domainSet.begin();
          const std::set<size_t>::const_iterator lastDomain = domainSet.end();

          while (domainItr != lastDomain)
            {
              partsItr = firstPart;
              while (partsItr != lastPart)
                {
                  domainGhostIndex[*domainItr].push_back(*partsItr);
                  partsItr++;
                }
              domainItr++;
            }
        }
    }

  return domainGhostIndex;
}

///
/// calculates the center of the universe
/// yes, it can be done :-)
///
void CostZone::centerOfTheUniverse(void)
{
  matrixRefType pos(PartManager.pos);
  bitsetRefType blk(PartManager.blacklisted);
  const size_t noParts = PartManager.getNoLocalParts();

  valueType xMin = std::numeric_limits<valueType>::max();
  valueType yMin = std::numeric_limits<valueType>::max();
  valueType zMin = std::numeric_limits<valueType>::max();

  valueType xMax = std::numeric_limits<valueType>::min();
  valueType yMax = std::numeric_limits<valueType>::min();
  valueType zMax = std::numeric_limits<valueType>::min();

  for (size_t i = 0; i < noParts; i++)
    {
      ///
      /// if particle is blacklisted, skip it
      ///
      if (blk[i])
        continue;

      xMin = pos(i, X) < xMin ? pos(i, X) : xMin;
      xMax = pos(i, X) > xMax ? pos(i, X) : xMax;

      yMin = pos(i, Y) < yMin ? pos(i, Y) : yMin;
      yMax = pos(i, Y) > yMax ? pos(i, Y) : yMax;

      zMin = pos(i, Z) < zMin ? pos(i, Z) : zMin;
      zMax = pos(i, Z) > zMax ? pos(i, Z) : zMax;
    }

#ifdef SPHLATCH_PARALLEL
  CommManager.min(xMin);
  CommManager.max(xMax);

  CommManager.min(yMin);
  CommManager.max(yMax);

  CommManager.min(zMin);
  CommManager.max(zMax);
#endif

  ///
  /// the box size is increased by a small factor, so that the cell
  /// indices in the calculation below never become negative or >= noCells1D
  /// due to rounding errors
  ///
#ifdef SPHLATCH_SINGLEPREC
  sidelength = (1 + 1e-4) *
#else
  sidelength = (1 + 1e-8) *
#endif
               std::max(std::max(xMax - xMin, yMax - yMin), zMax - zMin);

  xCenter = (xMin + xMax) / 2.;
  yCenter = (yMin + yMax) / 2.;
  zCenter = (zMin + zMax) / 2.;
}

///
/// put the local particles in the costzone cells
///
void CostZone::fillCostzoneCells()
{
  matrixRefType pos(PartManager.pos);
  bitsetRefType blk(PartManager.blacklisted);
  const size_t noParts = PartManager.getNoLocalParts();

  ///
  /// empty the costzone cells
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
  const valueType lengthToIndex = static_cast<valueType>(noCells1D) / sidelength;

  ///
  /// put the local particles in the corresponding costzone cells
  ///
  for (size_t i = 0; i < noParts; i++)
    {
      ///
      /// if particle is blacklisted, skip it
      ///
      if (blk[i])
        continue;

      const size_t xIndex = lrint((pos(i, X) - xMin) * lengthToIndex - 0.5);
      const size_t yIndex = lrint((pos(i, Y) - yMin) * lengthToIndex - 0.5);
      const size_t zIndex = lrint((pos(i, Z) - zMin) * lengthToIndex - 0.5);

      ///
      /// if this assertion fails, increase the factor in front
      /// of the sidelength calculation above
      ///
      assert(xIndex < noCells1D);
      assert(yIndex < noCells1D);
      assert(zIndex < noCells1D);

      const size_t cartIndex = xIndex + yIndex * noCells1D + zIndex * noCells2D;

      costzoneCells[cartIndex].push_back(i);
      partsPerCell[cartIndex] += 1;
    }
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

valueType CostZone::getAtomicLength(void)
{
  return (sidelength / static_cast<valueType>(noCells1D));
}

size_t CostZone::getDepth(void)
{
  return depth;
}

size_t CostZone::getNoGhosts(void)
{
  return noGhosts;
}
};

#endif
