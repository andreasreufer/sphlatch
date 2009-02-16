#ifndef BHTREE_DYNAMIC_H
#define BHTREE_DYNAMIC_H

/*
 *  bhtree_dynamic.h
 *
 *  base header file for the dynamic SPHLATCH Barnes&Hut tree
 *
 *  Created by Andreas Reufer on 01.12.08.
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#include <omp.h>

#include "typedefs.h"

#include "bhtree_node_cells.h"
#include "bhtree_node_particle.h"

#include "allocator_simple.h"
#include "allocator_chunk.h"

#include "particle_manager.h"
#include "communication_manager.h"

namespace sphlatch {
class BHTree {
public:
   friend class BHTreeWorker;
   friend class BHTreeCZBuilder;
   friend class BHTreePartsInsertMover;
   friend class BHTreeWorkerGrav;

   class BHTreeConfig {
public:
      BHTreeConfig() :
         partAllocSize(1048576),
         cellAllocSize(131072),
         czllAllocSize(1024),
         maxNoParts(1048576),
         overSize(1.5)
      { }

      ~BHTreeConfig() { }
      size_t partAllocSize, cellAllocSize, czllAllocSize;
      size_t maxNoParts;

      fType xMin, yMin, zMin;
      fType xMax, yMax, zMax;

      fType overSize;
   };

   typedef BHTree                         selfType;
   typedef BHTree&                        selfRef;
   typedef BHTree*                        selfPtr;

   typedef genericNode*                   nodePtrT;
   typedef const genericNode*             nodePtrCT;

   typedef particleNode                   partT;
   typedef particleNode*                  partPtrT;

   typedef genericCellNode                gcllT;
   typedef genericCellNode*               gcllPtrT;

   typedef quadrupoleCellNode             cellT;
   typedef quadrupoleCellNode*            cellPtrT;

   typedef costzoneCellNode               czllT;
   typedef costzoneCellNode*              czllPtrT;

   typedef std::list<czllPtrT>            czllPtrListT;
   typedef czllPtrListT::iterator         czllPtrListItrT;
   typedef czllPtrListT::const_iterator   czllPtrListCItrT;

   typedef std::vector<partPtrT>          partPtrVectT;

   BHTree(BHTreeConfig _config);
   ~BHTree();
   static selfRef instance(BHTreeConfig _config);

private:
   static selfPtr _instance;

protected:
   nodePtrT rootPtr;

   ChunkAllocator<partT>         partAllocator;
   ChunkAllocator<cellT>         cellAllocator;
   ChunkAllocator<czllT>         czllAllocator;

   czllPtrListT CZbottomCells;

   partPtrVectT partProxies;

   const BHTreeConfig config;

private:
   typedef ParticleManager        partManagerT;
   typedef CommunicationManager   commManagerT;

   partManagerT& PartManager;
   commManagerT& CommManager;

   matrixRefType pos;
};

BHTree::BHTree(BHTreeConfig _config) :
   partAllocator(_config.partAllocSize),
   cellAllocator(_config.cellAllocSize),
   czllAllocator(_config.czllAllocSize),
   partProxies(_config.maxNoParts),
   config(_config),
   PartManager(partManagerT::instance()),
   CommManager(commManagerT::instance()),
   pos(PartManager.pos)
{
   ///
   /// allocate root cell and set cell
   /// size, add root cell to CZbottom
   /// cells list
   ///
   rootPtr = czllAllocator.pop();
   static_cast<czllPtrT>(rootPtr)->clear();

   ///
   /// this cell is a bottom cell, add it to the list
   ///
   CZbottomCells.push_back(static_cast<czllPtrT>(rootPtr));
   static_cast<czllPtrT>(rootPtr)->listItr  = CZbottomCells.end();
   static_cast<czllPtrT>(rootPtr)->atBottom = true;

   ///
   /// determine the root cell size
   ///
   fType xMin = std::numeric_limits<fType>::max();
   fType yMin = std::numeric_limits<fType>::max();
   fType zMin = std::numeric_limits<fType>::max();

   fType xMax = std::numeric_limits<fType>::min();
   fType yMax = std::numeric_limits<fType>::min();
   fType zMax = std::numeric_limits<fType>::min();

   using namespace vectindices;

   const size_t noParts = PartManager.getNoLocalParts();
   for (size_t i = 0; i < noParts; i++)
   {
      xMin = pos(i, X) < xMin ? pos(i, X) : xMin;
      xMax = pos(i, X) > xMax ? pos(i, X) : xMax;

      yMin = pos(i, Y) < yMin ? pos(i, Y) : yMin;
      yMax = pos(i, Y) > yMax ? pos(i, Y) : yMax;

      zMin = pos(i, Z) < zMin ? pos(i, Z) : zMin;
      zMax = pos(i, Z) > zMax ? pos(i, Z) : zMax;
   }

   CommManager.min(xMin);
   CommManager.max(xMax);
   CommManager.min(yMin);
   CommManager.max(yMax);
   CommManager.min(zMin);
   CommManager.max(zMax);

   static_cast<czllPtrT>(rootPtr)->xCen  = 0.5*(xMin + xMax);
   static_cast<czllPtrT>(rootPtr)->yCen  = 0.5*(yMin + yMax);
   static_cast<czllPtrT>(rootPtr)->zCen  = 0.5*(zMin + zMax);
   static_cast<czllPtrT>(rootPtr)->clSz  = config.overSize *
                         std::max(std::max(xMax - xMin,
                                           yMax - yMin), zMax - zMin);

   ///
   /// resize particle proxy vector
   ///
   const size_t maxNoParts = config.maxNoParts;
   for (size_t i = 0; i < maxNoParts; i++)
   {
      partProxies[i] = NULL;
   }
}

BHTree::~BHTree()
{ }

BHTree::selfPtr BHTree::_instance = NULL;
BHTree::selfRef BHTree::instance(BHTreeConfig _config)
{
   if (_instance == NULL)
      _instance = new BHTree(_config);
   return(*_instance);
}
};
#endif
