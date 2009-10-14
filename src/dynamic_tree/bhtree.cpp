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

#include "bhtree_nodes.cpp"
#include "bhtree.h"

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
   
   btStart = static_cast<czllPtrT>(rootPtr);
   btEnd   = static_cast<czllPtrT>(rootPtr);
   lbtStart = NULL;
   lbtEnd   = NULL;

   static_cast<czllPtrT>(rootPtr)->clear();
   static_cast<czllPtrT>(rootPtr)->atBottom = true;
   static_cast<czllPtrT>(rootPtr)->depth    = 0;
   static_cast<czllPtrT>(rootPtr)->ident    = 0;
   noCells++;

   // just temporary
   static_cast<czllPtrT>(rootPtr)->cen  = 0.5, 0.5, 0.5;
   static_cast<czllPtrT>(rootPtr)->clSz = 1.;

   maxDepth = 100;

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

};
#endif
