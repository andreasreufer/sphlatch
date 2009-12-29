#ifndef BHTREE_H
#define BHTREE_H

/*
 *  bhtree.h
 *
 *  header file for the dynamic SPHLATCH Barnes&Hut tree
 *
 *  Created by Andreas Reufer on 07.10.09.
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"

#include "bhtree_nodes.h"
#include "bhtree_housekeeper.h"
#include "bhtree_part_insertmover.h"

namespace sphlatch {

class BHTree {
public:
   friend class BHTreeWorker;
   friend class BHTreeCZBuilder;
   friend class BHTreePartsInsertMover;
   friend class BHTreeWorkerGrav;

   typedef BHTree                   selfType;
   typedef BHTree&                  selfRef;
   typedef BHTree*                  selfPtr;

   typedef std::list<czllPtrT>      czllPtrListT;
   typedef std::vector<czllPtrT>    czllPtrVectT;

   BHTree();
   ~BHTree();

   static selfRef instance(void);

   ///
   /// some constants
   ///
   static const size_t maxDepth       = 128;
   static const size_t maxCZBottCells = 16384;

   static const fType cellsPerThread = 16;

   ///
   /// public functions
   ///
   void setExtent(const box3dT _box);
   void insertPart(treeGhost& _part);
   void update(const fType _cmin, const fType _cmax);
   czllPtrVectT getCZbottomLoc();

private:
   static selfPtr _instance;

   czllPtrVectT getCzllPtrVect(czllPtrListT _czllList);
   void sumUpCosts(), sumUpCostsRec();
   size_t round;

protected:
   nodePtrT rootPtr;

   ///
   /// first and last (local) CZ cell at bottom
   ///
   std::list<czllPtrT> CZbottom, CZbottomLoc;

   size_t noCells, noParts;

   const size_t noThreads;

private:
   BHTreePartsInsertMover insertmover;
};
};
#endif
