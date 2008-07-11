#ifndef SPHLATCH_BHTREE_ERRHANDLER_H
#define SPHLATCH_BHTREE_ERRHANDLER_H

/*
 *  bhtree_errhandler.h
 *
 *
 *  Created by Andreas Reufer on 10.07.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"
#include "err_handler.h"

#include "bhtree_generic_node.h"
#include "bhtree_cell_node.h"


namespace sphlatch
{
class TreeTooDeep : public GenericError
{
public:
typedef genericNode* nodePtrT;
typedef genericCellNode* cellPtrT;

TreeTooDeep(size_t _depth, size_t _part, nodePtrT _rootPtr);
~TreeTooDeep();
};

TreeTooDeep::TreeTooDeep(size_t _depth, size_t _part, nodePtrT _rootPtr)
{
  using namespace sphlatch::vectindices;
  matrixRefType pos(PartManager.pos);

  Logger.stream << "error: tree too deep, part " << _part << " ["
                << pos(_part, X) << ","
                << pos(_part, Y) << ","
                << pos(_part, Z) << "]"
                << " insert at depth " << _depth
                << " (tree root: ["
                << static_cast<cellPtrT>(_rootPtr)->xCenter << ","
                << static_cast<cellPtrT>(_rootPtr)->yCenter << ","
                << static_cast<cellPtrT>(_rootPtr)->zCenter << "] size "
                << static_cast<cellPtrT>(_rootPtr)->cellSize << ")";
  Logger.flushStream();
  Logger.destroy();
};

TreeTooDeep::~TreeTooDeep()
{
};
};

#endif
