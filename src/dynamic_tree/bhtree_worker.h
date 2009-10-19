#ifndef BHTREE_WORKER_H
#define BHTREE_WORKER_H

/*
 *  bhtree_worker.h
 *
 *  Created by Andreas Reufer on 07.10.09.
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#include "bhtree_nodes.h"

#include "bhtree.h"
//class BHTree;

namespace sphlatch {
class BHTreeWorker {
public:
   typedef BHTree&                    treeRefT;
   typedef BHTree*                    treePtrT;

   //typedef BHTree::czllPtrListT    czllPtrListT;
   typedef std::list<czllPtrT> czllPtrListT;

   BHTreeWorker(const BHTreeWorker& _worker);
   BHTreeWorker(const treePtrT _treePtr);
   ~BHTreeWorker();

protected:
   void goRoot();
   void goUp();
   void goNext();
   void goSkip();
   void goChild(const size_t _n);
   void goNeighbour(const size_t _n);

   bool pointInsideCell(const vect3dT& _pos);
   bool pointInsideCell(const vect3dT& _pos, const nodePtrT _nodePtr);

   size_t getOctant(const vect3dT& _pos);
   size_t getOctant(const vect3dT& _pos, const nodePtrT _nodePtr);

   size_t getChildNo(nodePtrT _nodePtr);

   void czllToCell(nodePtrT _nodePtr);
   void cellToCZll(nodePtrT _nodePtr);
   void partToCell(nodePtrT _nodePtr, const size_t _oct);

   const size_t   noThreads, myThread;
   const treePtrT treePtr;

   nodePtrT curPtr, rootPtr;
};


///
/// specialization for a read-only BHTree worker (eg. gravity walk)
///
class BHTreeWorkerRO : public BHTreeWorker {
public:
   BHTreeWorkerRO(treePtrT _treePtr) : BHTreeWorker(_treePtr) { }
   BHTreeWorkerRO(const BHTreeWorkerRO& _workerRO) : BHTreeWorker(_workerRO) { }

protected:
   ///
   /// overwrite curPtr by a faster read-only version
   ///
   nodePtrT curPtr, rootPtr;
};
};
#endif
