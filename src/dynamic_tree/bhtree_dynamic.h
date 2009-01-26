#ifndef BHTREE_DYNAMIC_H
#define BHTREE_DYNAMIC_H

/*
 *  bhtree_dynamic.h
 *
 *  base header file for the dynamic SPHLATCH Barnes&Hut tree
 *
 *  Created by Andreas Reufer on 01.12.08.
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#include <omp.h>

#include "typedefs.h"

#include "bhtree_node_cells.h"
#include "bhtree_node_particle.h"
#include "stack_allocator.h"

namespace sphlatch {
class BHTree {
public:
   friend class BHTreeWorker;
   friend class BHTreeCZBuilder;
   friend class BHTreePartsInsertMover;

   typedef BHTree                  selfType;
   typedef BHTree&                 selfRef;
   typedef BHTree*                 selfPtr;

   typedef genericNode*            nodePtrT;

   typedef particleNode            partT;
   typedef particleNode*           partPtrT;

   typedef quadrupoleCellNode      cellT;
   typedef quadrupoleCellNode*     cellPtrT;

   typedef costzoneCellNode        czllT;
   typedef costzoneCellNode*       czllPtrT;

   typedef std::list<czllPtrT>     czllPtrListT;

   typedef std::vector<partPtrT>   partPtrVectT;

   BHTree();
   ~BHTree();

   static selfRef instance(void);

private:
   static selfPtr _instance;

protected:
   czllPtrT rootPtr;

   StackAllocator<partT, 1048576> partAllocator;
   StackAllocator<cellT, 262144>  cellAllocator;
   StackAllocator<czllT, 16384>   czllAllocator;

   czllPtrListT CZbottomCells;

   partPtrVectT partProxies;
};

BHTree::BHTree()
{
   ///
   /// allocate root cell and set cell
   /// size, add root cell to CZbottom
   /// cells list
   ///
   rootPtr           = czllAllocator.pop();
   rootPtr->atBottom = true;

   ///
   /// resize particle proxy vector
   ///
   partProxies.reserve(1048576);
   
   partProxies.resize(1048576);
   for (size_t i = 0; i < 1048576; i++)
   {
     partProxies[i] = NULL;
   }
   partProxies.resize(0);
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
};
#endif
