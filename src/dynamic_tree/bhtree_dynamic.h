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

#ifdef SPHLATCH_OPENMP
 #include <omp.h>
#endif

#include "typedefs.h"


#include "bhtree_node_cells.h"
#include "bhtree_node_particle.h"
#include "bhtree_particle.h"

namespace sphlatch {
class BHTree {
public:
   friend class BHTreeWorker;
   friend class BHTreeCZBuilder;
   friend class BHTreePartsInsertMover;
   friend class BHTreeWorkerGrav;

   typedef BHTree                selfType;
   typedef BHTree&               selfRef;
   typedef BHTree*               selfPtr;

   typedef genericNode*          nodePtrT;
   typedef const genericNode*    nodePtrCT;

   typedef treeGhost             partT;
   typedef treeGhost*            partPtrT;

   typedef particleNode          pnodT;
   typedef particleNode*         pnodPtrT;

   typedef genericCellNode       gcllT;
   typedef genericCellNode*      gcllPtrT;

   typedef quadrupoleCellNode    cellT;
   typedef quadrupoleCellNode*   cellPtrT;

   typedef costzoneCellNode      czllT;
   typedef costzoneCellNode*     czllPtrT;

   BHTree();
   ~BHTree();

   static selfRef instance(void);

private:
   static selfPtr _instance;

protected:
   nodePtrT rootPtr;

   ///
   /// first and last (local) CZ cell at bottom
   ///
   czllPtrT btStart, btEnd;
   czllPtrT lbtStart, lbtEnd;

   size_t noCells, noParts;

   const size_t noThreads;
};

BHTree::BHTree() :
   noCells(0),
   noParts(0),
   noThreads(omp_get_num_threads())
{
   ///
   /// allocate root cell and set cell
   /// size, add root cell to CZbottom
   /// cells list
   ///
   rootPtr = new czllT;

   btStart  = static_cast<czllPtrT>(rootPtr);
   btEnd    = static_cast<czllPtrT>(rootPtr);
   lbtStart = NULL;
   lbtEnd   = NULL;

   static_cast<czllPtrT>(rootPtr)->clear();
   static_cast<czllPtrT>(rootPtr)->atBottom = true;
   static_cast<czllPtrT>(rootPtr)->depth    = 0;
   static_cast<czllPtrT>(rootPtr)->ident    = 0;
   noCells++;

   //FIXME: just temporary
   static_cast<czllPtrT>(rootPtr)->cen  = 0.5, 0.5, 0.5;
   static_cast<czllPtrT>(rootPtr)->clSz = 1.;
}

BHTree::~BHTree()
{ }

BHTree::selfPtr BHTree::_instance = NULL;
BHTree::selfRef BHTree::instance()
{
   if (_instance == NULL)
      _instance = new BHTree;
   return(*_instance);
}
};
#endif
