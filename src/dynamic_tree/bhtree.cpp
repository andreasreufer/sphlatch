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
#include "bhtree_cz_builder.cpp"

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

   CZbottom.push_back(static_cast<czllPtrT>(rootPtr));

   static_cast<czllPtrT>(rootPtr)->clear();
   static_cast<czllPtrT>(rootPtr)->atBottom = true;
   static_cast<czllPtrT>(rootPtr)->depth    = 0;
   static_cast<czllPtrT>(rootPtr)->ident    = 0;
   noCells++;

   // just temporary
   static_cast<czllPtrT>(rootPtr)->cen  = 0.5, 0.5, 0.5;
   static_cast<czllPtrT>(rootPtr)->clSz = 1.;

   static_cast<czllPtrT>(rootPtr)->parent = NULL;

   //maxDepth = 100;

/*   std::cout << static_cast<czllPtrT>(rootPtr)->cen << "  "
             << static_cast<czllPtrT>(rootPtr)->clSz << "\n";*/
}

BHTree::~BHTree()
{ }

BHTree::selfPtr BHTree::_instance = NULL;
BHTree::selfRef BHTree::instance()
{
   if (_instance == NULL)
      _instance = new BHTree;
   return(*_instance);
}

void BHTree::insertParts(partVectT& _parts)
{
   const size_t noParts = _parts.size();

   BHTreePartsInsertMover inserter(this);

   for (size_t i = 0; i < noParts; i++)
   {
      inserter.insert(_parts[i]);
   }
}

void BHTree::update()
{
   // move particles
   // (prepare next walk?)
   BHTreePartsInsertMover       mover(this);
   czllPtrListT::iterator       CZItr = CZbottom.begin();
   czllPtrListT::const_iterator CZEnd = CZbottom.end();

   while (CZItr != CZEnd )
   {
     std::cout << *CZItr << "\n";
     mover.move(*CZItr);
     CZItr++;
   } 

   // rebalance trees
   BHTreeCZBuilder czbuilder(this);
   czbuilder.rebalance();

   // exchange costzone cells and their particles

   // push down orphans
   while (CZItr != CZEnd )
   {
     std::cout << *CZItr << "\n";
     mover.move(*CZItr);
     CZItr++;
   }

   // clean up tree

   // prepare walks (next & skip)

   std::cout << "prepare walks \n";
   omp_set_num_threads(4);

   czllPtrVectT CZBottomV = getCzllPtrVect(CZbottom);
   const int noCZBottomCells = CZBottomV.size();

   BHTreeHousekeeper HK(this);
#pragma omp parallel for firstprivate(HK)
   for (int i = 0; i < noCZBottomCells; i++)
   {
     HK.setNext( CZBottomV[i] );
     std::cout << i << " " << &HK << " " << omp_get_thread_num() << "\n";
   }
   HK.setNextCZ();
   HK.setSkip();

   // exchange MP moments
}

BHTree::czllPtrVectT BHTree::getCzllPtrVect(czllPtrListT _czllList)
{
  czllPtrVectT vect;
  vect.resize(_czllList.size());

  czllPtrListT::iterator       CZItr = _czllList.begin();
  czllPtrListT::const_iterator CZEnd = _czllList.end();

  size_t i = 0;
  while (CZItr != CZEnd)
  {
    vect[i] = *CZItr;
    CZItr++;
    i++;
  }
  return vect;
}
};
#endif
