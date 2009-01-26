#ifndef BHTREE_PART_INSERTMOVER_H
#define BHTREE_PART_INSERTMOVER_H

/*
 *  bhtree_part_insertmover.h
 *
 *  Created by Andreas Reufer on 02.12.08.
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#include "bhtree_worker.h"
#include "bhtree_errhandler_tmp.h"

namespace sphlatch {
class BHTreePartsInsertMover : public BHTreeWorker {
public:

   typedef BHTreePartsInsertMover   selfT;

   BHTreePartsInsertMover(treePtrT _treePtr) 
     : BHTreeWorker(_treePtr) { }
   BHTreePartsInsertMover(const selfT& _inserter) 
     : BHTreeWorker(_inserter) { }
   ~BHTreePartsInsertMover() { }

public:
   void insert(const size_t _i);
   void moveAll();

private:
   void recursor(const size_t _i);
   void move(const size_t _i);
};

///
/// entry function to insert particle:
/// - check whether particle lies inside root cell
/// - if so, call insertion recursor
///
void BHTreePartsInsertMover::insert(const size_t _i)
{
   goRoot();

   ///
   /// allocate particle node, wire it to its proxy and
   /// set nodes parameters
   ///
   const partPtrT curPartPtr = treePtr->partAllocator.pop();
   treePtr->partProxies[_i] = curPartPtr;

   static_cast<partPtrT>(curPartPtr)->ident = _i;
   static_cast<partPtrT>(curPartPtr)->depth = 0;

   ///
   /// wire the root cell as particles parent (in case the tree
   /// is not cleaned before use, this particle is not found in
   /// a tree walk, as it is not yet wired as a child of a cell)
   ///
   static_cast<partPtrT>(curPartPtr)->parent = treePtr->rootPtr;
}


};

#endif
