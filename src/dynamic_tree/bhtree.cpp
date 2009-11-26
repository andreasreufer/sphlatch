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
#include "bhtree_worker_mp.cpp"

#include "bhtree_treedump.cpp"

namespace sphlatch {
BHTree::BHTree() :
   rootPtr(new czllT),
   noCells(1),
   noParts(0),
#ifdef SPHLATCH_OPENMP
   noThreads(omp_get_num_threads()),
#else
   noThreads(1),
#endif
   insertmover(this)
{
   ///
   /// allocate root cell and set cell
   /// size, add root cell to CZbottom
   /// cells list
   ///
   CZbottom.push_back(static_cast<czllPtrT>(rootPtr));

   static_cast<czllPtrT>(rootPtr)->clear();
   static_cast<czllPtrT>(rootPtr)->atBottom = true;
   static_cast<czllPtrT>(rootPtr)->depth    = 0;
   static_cast<czllPtrT>(rootPtr)->ident    = 0;

   // just temporary
   static_cast<czllPtrT>(rootPtr)->cen  = 0.5, 0.5, 0.5;
   static_cast<czllPtrT>(rootPtr)->clSz = 1.;

   static_cast<czllPtrT>(rootPtr)->parent = NULL;

   round = 0;

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

void BHTree::insertPart(treeGhost& _part)
{
   insertmover.insert(_part);
}

void BHTree::update(const fType _cmarkLow, const fType _cmarkHigh)
{
   //std::cout << "entered Tree.update() ... round " << round << "\n";

   BHTreeDump         dumper(this);
   std::ostringstream roundStr;
   roundStr << round;
   std::string dumpName = "dump";
   dumpName.append(roundStr.str());

   // rebalance trees
   //dumper.dotDump(dumpName + "_0.dot");
   //dumper.ptrDump(dumpName + "_0.ptr");

   // move particles
   // (prepare next walk?)
   czllPtrListT::iterator       CZItr = CZbottom.begin();
   czllPtrListT::const_iterator CZEnd = CZbottom.end();
   std::cout << "move parts in " << CZbottom.size() << " CZ cells\n";
   while (CZItr != CZEnd)
   {
      insertmover.move(*CZItr);
      CZItr++;
   }

   // rebalance trees
   //dumper.dotDump(dumpName + "_1.dot");
   //dumper.ptrDump(dumpName + "_1.ptr");

   const fType normCellCost = 1. / ( noThreads * cellsPerThread );
   const fType costMin = normCellCost * _cmarkLow;
   const fType costMax = normCellCost * _cmarkHigh;

   BHTreeCZBuilder czbuilder(this);
   czbuilder.rebalance(costMin, costMax);

   //dumper.dotDump(dumpName + "_2.dot");
   //dumper.ptrDump(dumpName + "_2.txt");

   // compose vector of CZ cell pointers
   czllPtrVectT CZbottomV       = getCzllPtrVect(CZbottom);
   const int    noCZBottomCells = CZbottomV.size();

   // exchange costzone cells and their particles

   // push down orphans
   //std::cout << "\npush down orphans\n";
   CZItr = CZbottom.begin();
   while (CZItr != CZEnd)
   {
      insertmover.pushDownOrphans(*CZItr);
      CZItr++;
   }

   //dumper.dotDump(dumpName + "_3.dot");
   //dumper.ptrDump(dumpName + "_3.txt");

   // clean up tree

   // prepare walks (next & skip)
   //std::cout << "housekeeping          \n";

   BHTreeHousekeeper HK(this);
   BHTreeMPWorker    MP(this);
#pragma omp parallel for firstprivate(HK, MP)
   for (int i = 0; i < noCZBottomCells; i++)
   {
      // set next pointers
      HK.setNext(CZbottomV[i]);

      // clean up
      HK.minTree(CZbottomV[i]);

      // calculate MP moments
      MP.calcMultipoles(CZbottomV[i]);
   }
   //dumper.dotDump(dumpName + "_4.dot");
   //dumper.ptrDump(dumpName + "_4.txt");

   //std::cout << "prepare CZ next walk  \n";
   HK.setNextCZ();
   //std::cout << "prepare    skip walk  \n";
   HK.setSkip();
   //std::cout << "calculate MP moments  \n";
   MP.calcMultipolesCZ();

   //dumper.dotDump(dumpName + "_5.dot");
   //dumper.ptrDump(dumpName + "_5.txt");

   // exchange MP moments
   round++;

   //FIXME: change this in parallel version
   CZbottomLoc = CZbottom;
}

BHTree::czllPtrVectT BHTree::getCZbottomLoc()
{
   return(getCzllPtrVect(CZbottomLoc));
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
   return(vect);
}
};
#endif
