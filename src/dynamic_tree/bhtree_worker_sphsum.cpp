#ifndef BHTREE_WORKER_SPHSUM_CPP
#define BHTREE_WORKER_SPHSUM_CPP

/*
 *  bhtree_worker_sphsum.cpp
 *
 *  Created by Andreas Reufer on 14.10.09.
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#include "bhtree_worker.cpp"
#include "bhtree_particle.h"

namespace sphlatch {
template<typename Summand>
class SPHsumWorker : public BHTreeWorker {
public:

   SPHsumWorker(const treePtrT _treePtr) : BHTreeWorker(_treePtr) { }
   SPHsumWorker(const BHTree& _) : BHTreeWorker(_) { }
   ~SPHsumWorker() { }

   Summand S;

   void operator()()
   {
      treeGhost p1, p2;

      S(&p1, &p2);
   }
};
};

#endif
