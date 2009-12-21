#ifndef RANKLIST_H
#define RANKLIST_H

/*
 *  ranklist.h
 *
 *
 *  Created by Andreas Reufer on 10.04.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"

#include "particle.h"
#include "memorymanager.h"

namespace sphlatch {
class Ranklist {
/*typedef genericNode nodeT;

   typedef genericNode* nodePtrT;
   typedef particleNode* partPtrT;
   typedef genericCellNode* cellPtrT;
   typedef particleProxy* partProxyPtrT;*/

typedef sphlatch::MemoryManager memManagerType;
memManagerType& MemManager;
matrixRefType Data;

///
/// constructor
///
public:
Ranklist() :
  MemManager(memManagerType::instance()),
  Data(MemManager.Data)
{
  partsPerCell = 5;
  neighbourList.resize(10000);
};

///
/// destructor
///
~Ranklist(void)
{
}

private:
partsIndexVectType indicesRankedX, indicesRankedY, indicesRankedZ,
                   rsIndicesX, rsIndicesY, rsIndicesZ;

// replace later on by valVectType
std::valarray<fType> cellMinX, cellMinY, cellMinZ,
                         cellMaxX, cellMaxY, cellMaxZ;

size_t noParts, noCells, partsPerCell;

public:
//void insertParticle(partProxyPtrT _newPayload)
//   {
//   }

///
/// small helper class to sort the rankspace index vectors
///
private:
template <size_t _idx> class smallerThan {
typedef sphlatch::MemoryManager memManagerType;
memManagerType& MemManager;

public:
smallerThan(void) :
  MemManager(memManagerType::instance())
{
}

bool operator()(const size_t& _i, const size_t& _j)
{
  return(MemManager.Data(_i, _idx) < MemManager.Data(_j, _idx));
}
};

public:
void prepare()
{
  noParts = Data.size1();

  indicesRankedX.resize(noParts);
  indicesRankedY.resize(noParts);
  indicesRankedZ.resize(noParts);

  rsIndicesX.resize(noParts);
  rsIndicesY.resize(noParts);
  rsIndicesZ.resize(noParts);

  for (size_t i = 0; i < noParts; i++)
    {
      indicesRankedX[i] = i;
      indicesRankedY[i] = i;
      indicesRankedZ[i] = i;
    }

  std::sort(indicesRankedX.begin(), indicesRankedX.end(), smallerThan<X>());
  std::sort(indicesRankedY.begin(), indicesRankedY.end(), smallerThan<Y>());
  std::sort(indicesRankedZ.begin(), indicesRankedZ.end(), smallerThan<Z>());

  ///
  /// calculate the inverse of the ranked indices
  ///
  for (size_t i = 0; i < noParts; i++)
    {
      rsIndicesX[ indicesRankedX[i] ] = i;
      rsIndicesY[ indicesRankedY[i] ] = i;
      rsIndicesZ[ indicesRankedZ[i] ] = i;
    }


  ///
  /// prepare the cell min/max containers
  ///
  noCells = lrint(ceil(static_cast<double>(noParts) / static_cast<double>(partsPerCell)));
  cellMinX.resize(noCells);
  cellMinY.resize(noCells);
  cellMinZ.resize(noCells);
  cellMaxX.resize(noCells);
  cellMaxY.resize(noCells);
  cellMaxZ.resize(noCells);

  ///
  /// determine the cell min/max
  ///
  size_t minRank, maxRank;
  for (size_t i = 0; i < noCells; i++)
    {
      minRank = i * partsPerCell;
      maxRank = std::min(minRank + partsPerCell - 1, noParts - 1);

      cellMinX[i] = Data(indicesRankedX[minRank], X);
      cellMaxX[i] = Data(indicesRankedX[maxRank], X);
      cellMinY[i] = Data(indicesRankedY[minRank], Y);
      cellMaxY[i] = Data(indicesRankedY[maxRank], Y);
      cellMinZ[i] = Data(indicesRankedZ[minRank], Z);
      cellMaxZ[i] = Data(indicesRankedZ[maxRank], Z);
    }

  ///
  ///
  ///
  for (size_t i = 0; i < noParts; i++)
    {
      indicesRankedX[ rsIndicesX[i] ] = rsIndicesY[i];
      indicesRankedY[ rsIndicesY[i] ] = rsIndicesZ[i];
    }

  rsIndicesY.resize(0);
  rsIndicesZ.resize(0);
}

public:
partsIndexVectType neighbourList;
/*void findNeighbours(const partProxyPtrT _curParticle,
                    const fRefType _search_radius)
   {
   }*/

private:
fType curPartX, curPartY, curPartZ;
size_t minRankX, maxRankX, minRankY, maxRankY, minRankZ, maxRankZ;

public:
void findNeighbours(const size_t& _curPartIndex,
                    const fRefType _search_radius)
{
  curPartX = Data(_curPartIndex, X);
  curPartY = Data(_curPartIndex, Y);
  curPartZ = Data(_curPartIndex, Z);

  minRankX = rsIndicesX[_curPartIndex] - (rsIndicesX[_curPartIndex] % partsPerCell);
  maxRankX = minRankX + partsPerCell - 1;

  minRankY = rsIndicesY[_curPartIndex] - (rsIndicesY[_curPartIndex] % partsPerCell);
  maxRankY = minRankY + partsPerCell - 1;

  minRankZ = rsIndicesZ[_curPartIndex] - (rsIndicesZ[_curPartIndex] % partsPerCell);
  maxRankZ = minRankZ + partsPerCell - 1;

  ///
  /// determine min. and max. rank in each dimension
  ///
  while ((curPartX - _search_radius) < (cellMinX[minRankX / partsPerCell]))
    {
      if (minRankX < partsPerCell)
        break;
      minRankX -= partsPerCell;
    }

  while ((curPartX + _search_radius) > (cellMaxX[((maxRankX + 1) / partsPerCell) - 1 ]))
    {
      if (maxRankX >= noParts - 1)
        break;
      maxRankX += partsPerCell;
    }

  while ((curPartY - _search_radius) < (cellMinY[minRankY / partsPerCell]))
    {
      if (minRankY < partsPerCell)
        break;
      minRankY -= partsPerCell;
    }

  while ((curPartY + _search_radius) > (cellMaxY[((maxRankY + 1) / partsPerCell) - 1 ]))
    {
      if (maxRankY >= noParts - 1)
        break;
      maxRankY += partsPerCell;
    }

  while ((curPartZ - _search_radius) < (cellMinZ[minRankZ / partsPerCell]))
    {
      if (minRankZ < partsPerCell)
        break;
      minRankZ -= partsPerCell;
    }

  while ((curPartZ + _search_radius) > (cellMaxZ[((maxRankZ + 1) / partsPerCell) - 1 ]))
    {
      if (maxRankZ >= noParts - 1)
        break;
      maxRankZ += partsPerCell;
    }

  static size_t curRankX, curRankY, curRankZ, neighIndex;
  static size_t noNeighbours;
  static fType neighDist;

  noNeighbours = 0;

  for (curRankX = minRankX; curRankX <= maxRankX; curRankX++)
    {
      curRankY = indicesRankedX[ curRankX ];
      if (curRankY < minRankY || curRankY > maxRankY)
        continue;

      curRankZ = indicesRankedY[ curRankY ];
      if (curRankZ < minRankZ || curRankZ > maxRankZ)
        continue;

      ///
      /// particle is in the same rankspace search volume like the particle,
      /// now check the distance between the two particles geometrically
      ///
      neighIndex = indicesRankedZ[ curRankZ ];

      neighDist = sqrt((Data(neighIndex, X) - curPartX)
                       * (Data(neighIndex, X) - curPartX)
                       + (Data(neighIndex, Y) - curPartY)
                       * (Data(neighIndex, Y) - curPartY)
                       + (Data(neighIndex, Z) - curPartZ)
                       * (Data(neighIndex, Z) - curPartZ));

      if (neighDist < _search_radius)
        {
          noNeighbours++;
          neighbourList[noNeighbours] = neighIndex;
        }
    }
  neighbourList[0] = noNeighbours;
}
};
};

#endif
