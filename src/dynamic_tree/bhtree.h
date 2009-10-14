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

   typedef BHTree                         selfType;
   typedef BHTree&                        selfRef;
   typedef BHTree*                        selfPtr;

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

   size_t maxDepth;



};

};
#endif
