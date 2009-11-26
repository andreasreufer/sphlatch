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


namespace sphlatch {
class BHTree;

class BHTreeWorker {
public:
   typedef BHTree&               treeRefT;
   typedef BHTree*               treePtrT;

   typedef std::list<czllPtrT>   czllPtrListT;

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
   
   bool sphereTotInCell(const vect3dT& _pos, const fType& _r);
   bool sphereTotOutCell(const vect3dT& _pos, const fType& _r);

   size_t getOctant(const vect3dT& _pos);
   size_t getOctant(const vect3dT& _pos, const nodePtrT _nodePtr);

   size_t getChildNo(nodePtrT _nodePtr);
   size_t getChildNo(nodePtrT _nodePtr, nodePtrT _parPtr);

   qcllPtrT czllToCell(nodePtrT _nodePtr);
   czllPtrT cellToCZll(nodePtrT _nodePtr);
   qcllPtrT partToCell(nodePtrT _nodePtr, const size_t _oct);

   const size_t   noThreads, myThread;
   const treePtrT treePtr;

   nodePtrT       curPtr;
   nodePtrT const rootPtr;
};


///
/// specialization for a read-only BHTree worker (eg. gravity walk)
///
//FIXME: doesn't work somehow
class BHTreeWorkerRO : public BHTreeWorker {
public:
   BHTreeWorkerRO(treePtrT _treePtr) : BHTreeWorker(_treePtr), curPtr(NULL) { }
   BHTreeWorkerRO(const BHTreeWorkerRO& _workerRO) : BHTreeWorker(_workerRO),
                                                     curPtr(NULL) { }

protected:
   ///
   /// overwrite curPtr by a faster (?) read-only version
   ///
   const nodePtrT curPtr;
};
};
#endif
