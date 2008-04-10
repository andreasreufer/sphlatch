#ifndef RANKSPACE_H
#define RANKSPACE_H

/*
 *  rankspace.h
 *
 *
 *  Created by Andreas Reufer on 10.04.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include "particle.h"
#include "typedefs.h"

//#include "communicationmanager.h"
#include "memorymanager.h"

namespace sphlatch {
class Rankspace {
/*typedef genericNode nodeT;

   typedef genericNode* nodePtrT;
   typedef particleNode* partPtrT;
   typedef genericCellNode* cellPtrT;
   typedef particleProxy* partProxyPtrT;*/

//typedef sphlatch::CommunicationManager commManagerType;
//commManagerType& CommManager;
typedef sphlatch::MemoryManager memManagerType;
memManagerType& MemManager;
matrixRefType Data;

///
/// constructor
///
public:
Rankspace() :
  MemManager(memManagerType::instance()),
  Data(MemManager.Data)
{
  neighbourList.resize(10000);
};

///
/// destructor
///
~Rankspace(void)
{
}

private:
partsIndexVectType indicesRankedX, indicesRankedY, indicesRankedZ,
                   rsIndicesX, rsIndicesY, rsIndicesZ;
size_t noParts;

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
    return ( MemManager.Data(_i, _idx) < MemManager.Data(_j, _idx) );
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

  std::sort( indicesRankedX.begin(), indicesRankedX.end(), smallerThan<X>() );
  std::sort( indicesRankedY.begin(), indicesRankedY.end(), smallerThan<Y>() );
  std::sort( indicesRankedZ.begin(), indicesRankedZ.end(), smallerThan<Z>() );

  for (size_t i = 0; i < noParts; i++)
  {
    rsIndicesX[ indicesRankedX[i] ] = i;
    rsIndicesY[ indicesRankedY[i] ] = i;
    rsIndicesZ[ indicesRankedZ[i] ] = i;
  }

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
                    const valueRefType _search_radius)
   {
   }*/

void findNeighbours(const size_t& _curPartIndex,
                    const valueRefType _search_radius)
{
}
};
};

#endif
