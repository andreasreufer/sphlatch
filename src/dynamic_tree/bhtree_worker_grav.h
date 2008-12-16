#ifndef BHTREE_WORKER_GRAV_H
#define BHTREE_WORKER_GRAV_H

/*
 *  bhtree_worker_grav.h
 *
 *  Created by Andreas Reufer on 14.12.08.
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#include "bhtree_worker.h"

namespace sphlatch {
class BHTreeWorkerGrav : public BHTreeWorker {
public:

/*typedef particleNode partType;
typedef quadrupoleCellNode cellType;
typedef costzoneCellNode czllType;

typedef genericNode* nodePtr;*/

///
/// constructors and destructors
///
BHTreeWorkerGrav(treePtrT _treePtr) : BHTreeWorker(_treePtr)
{
  const size_t tid = omp_get_thread_num();
  std::cout << tid << ":tpc:" << curPtr << ":" << this << "\n";
};

BHTreeWorkerGrav(const BHTreeWorkerGrav &_gravwork)
  : BHTreeWorker(_gravwork)
{
  const size_t tid = omp_get_thread_num();
  std::cout << tid << ":cpc:" << curPtr << ":" << this << "\n";
};
~BHTreeWorkerGrav() {};

};


};

#endif

