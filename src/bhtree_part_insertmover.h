#ifndef BHTREE_PART_INSERTMOVER_H
#define BHTREE_PART_INSERTMOVER_H

/*
 *  bhtree_part_insertmover.h
 *
 *  Created by Andreas Reufer on 01.11.09
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#include "bhtree_worker.h"
#include "bhtree_particle.h"

namespace sphlatch {
class BHTreePartsInsertMover : public BHTreeWorker {
public:
   typedef BHTreePartsInsertMover   selfT;
   typedef treeGhost                partT;

   BHTreePartsInsertMover(treePtrT _treePtr);
   BHTreePartsInsertMover(const selfT& _inserter);
   ~BHTreePartsInsertMover();

public:
   void insert(partT& _part);
   void move(const czllPtrT _czll);
   void pushDownOrphans(const czllPtrT _czll);

private:
   void pushUpAndToCZSingle(const pnodPtrT _pnodPtr);

   void pushDownSingle(const pnodPtrT _pnodPtr);

protected:
   const size_t maxDepth;
};
};

#endif
