#ifndef BHTREE_ERRHANDLER_H
#define BHTREE_ERRHANDLER_H

/*
 *  bhtree_errhandler.h
 *
 *  Created by Andreas Reufer on 16.12.08.
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"
#include "bhtree_node_cells.h"
#include "err_handler.h"

using namespace sphlatch::vectindices;

namespace sphlatch {
class PartOutsideTree : public GenericError
{
public:
   typedef genericNode*       nodePtrT;
   typedef genericCellNode*   cellPtrT;

   PartOutsideTree(size_t _newPart, nodePtrT _rootPtr);
   ~PartOutsideTree();
};


PartOutsideTree::PartOutsideTree(size_t   _newPart,
                                 nodePtrT _rootPtr)
{
   using namespace sphlatch::vectindices;
   matrixRefType pos(PartManager.pos);
   idvectRefType id(PartManager.id);
   const fType   rootSize = static_cast<cellPtrT>(_rootPtr)->clSz;
   const fType   rootX    = static_cast<cellPtrT>(_rootPtr)->xCen;
   const fType   rootY    = static_cast<cellPtrT>(_rootPtr)->yCen;
   const fType   rootZ    = static_cast<cellPtrT>(_rootPtr)->zCen;

#ifdef SPHLATCH_LOGGER
   Logger.stream
#else
   std::cerr
#endif
  << "error: partID "
  << id(_newPart) << " ["
  << pos(_newPart, X) << ","
  << pos(_newPart, Y) << ","
  << pos(_newPart, Z) << "]"
  << " outside tree (tree root: ["
  << rootX - 0.5 * rootSize << ","
  << rootY - 0.5 * rootSize << ","
  << rootZ - 0.5 * rootSize << "] - [ "
  << rootX + 0.5 * rootSize << ","
  << rootY + 0.5 * rootSize << ","
  << rootZ + 0.5 * rootSize << "] )";
#ifdef SPHLATCH_LOGGER
   Logger.flushStream();
   Logger.destroy();
#else
   std::cerr << "\n";
#endif
}

PartOutsideTree::~PartOutsideTree()
{ }
};

#endif
