#ifndef BHTREE_CPP
#define BHTREE_CPP

/*
 *  bhtree.cpp
 *
 *  Created by Andreas Reufer on 01.12.08.
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#ifdef SPHLATCH_OPENMP
 #include <omp.h>
#endif

#include "bhtree.h"
#include "bhtree_nodes.cpp"

#include "bhtree_housekeeper.cpp"
#include "bhtree_part_insertmover.cpp"

namespace sphlatch {

BHTree::BHTree() :
   noCells(0),
   noParts(0),
#ifdef SPHLATCH_OPENMP
   noThreads(omp_get_num_threads())
#else
   noThreads(1)
#endif
{
   ///
   /// allocate root cell and set cell
   /// size, add root cell to CZbottom
   /// cells list
   ///
   rootPtr = new czllT;
   
   CZbottom.push_back( static_cast<czllPtrT>(rootPtr) );

   static_cast<czllPtrT>(rootPtr)->clear();
   static_cast<czllPtrT>(rootPtr)->atBottom = true;
   static_cast<czllPtrT>(rootPtr)->depth    = 0;
   static_cast<czllPtrT>(rootPtr)->ident    = 0;
   noCells++;

   // just temporary
   static_cast<czllPtrT>(rootPtr)->cen  = 0.5, 0.5, 0.5;
   static_cast<czllPtrT>(rootPtr)->clSz = 1.;
   
   static_cast<czllPtrT>(rootPtr)->parent    = NULL;

   //maxDepth = 100;

/*   std::cout << static_cast<czllPtrT>(rootPtr)->cen << "  "
             << static_cast<czllPtrT>(rootPtr)->clSz << "\n";*/
};

BHTree::~BHTree()
{ };

BHTree::selfPtr BHTree::_instance = NULL;
BHTree::selfRef BHTree::instance()
{
   if (_instance == NULL)
      _instance = new BHTree;
   return(*_instance);
};

void BHTree::insertParts(partVectT& _parts)
{  
  const size_t noParts = _parts.size();

  BHTreePartsInsertMover inserter(this);

  for (size_t i = 0; i < noParts; i++)
  {
    inserter.insert(_parts[i]);
  }

};

void BHTree::update()
{
  // move particles
  // (prepare next walk?)
  BHTreePartsInsertMover mover(this);

  // rebalance trees
  // push down parts

  // exchange particles

  // push down orphans
  // clean up tree
  // prepare walks (next & skip)

  // exchange MP moments
};


};
#endif
