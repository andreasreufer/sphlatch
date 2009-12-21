#ifndef BHTREE_WORKER__H
#define BHTREE_WORKER__H

/*
 *  bhtree_worker_.h
 *
 *  Created by Andreas Reufer on 14.12.08.
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#include "bhtree_worker.h"

namespace sphlatch {
class BHTreeMyWorker : public BHTreeWorker {
public:

BHTreeMyWorker(const treePtrT _treePtr) : BHTreeWorker(_treePtr) {};
BHTreeMyWorker(const BHTree &_) : BHTreeWorker(_) {};
~BHTreeMyWorker() {};

};


};

#endif

