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

namespace sphlatch {
class BHTree {
public:
   friend class BHTreeWorker;
   friend class BHTreeCZBuilder;
   friend class BHTreePartsInsertMover;
   friend class BHTreeWorkerGrav;

   typedef BHTree    selfType;
   typedef BHTree&   selfRef;
   typedef BHTree*   selfPtr;

   typedef std::list<czllPtrT> czllPtrListT;

   BHTree();
   ~BHTree();

   static selfRef instance(void);

   ///
   /// some constants
   ///
   static const size_t maxDepth       = 128;
   static const size_t maxCZBottCells = 16384;

   ///
   /// public functions
   ///
   void insertParticles();
   void update();

private:
   static selfPtr _instance;

protected:
   nodePtrT rootPtr;

   ///
   /// first and last (local) CZ cell at bottom
   ///
   std::list<czllPtrT> CZbottom, CZbottomLoc;

   size_t noCells, noParts;

   const size_t noThreads;
};
};
#endif
