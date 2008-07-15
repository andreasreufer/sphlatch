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
class PartsTooClose : public GenericError
{
public:
typedef genericNode* nodePtrT;
typedef genericCellNode* cellPtrT;

PartsTooClose(size_t _depth, size_t _resPart, size_t _newPart, nodePtrT _rootPtr);
~PartsTooClose();
};

PartsTooClose::PartsTooClose(size_t _depth,
                                   size_t _resPart,
                                   size_t _newPart,
                                   nodePtrT _rootPtr)
{
  using namespace sphlatch::vectindices;
  matrixRefType pos(PartManager.pos);
  idvectRefType id(PartManager.id);
  const rootSize =  static_cast<cellPtrT>(_rootPtr)->cellSize;

  Logger.stream << "error: tree too deep (" << _depth << "), part ID "
                << id(_resPart) << " ["
                << pos(_resPart, X) << ","
                << pos(_resPart, Y) << ","
                << pos(_resPart, Z) << "]"
                << " and part ID "
                << id(_newPart) << " ["
                << pos(_newPart, X) << ","
                << pos(_newPart, Y) << ","
                << pos(_newPart, Z) << "]"
                << " too close "
                << " (tree root: ["
                << static_cast<cellPtrT>(_rootPtr)->xCenter - 0.5*rootSize
                << ","
                << static_cast<cellPtrT>(_rootPtr)->yCenter - 0.5*rootSize
                << ","
                << static_cast<cellPtrT>(_rootPtr)->zCenter - 0.5*rootSize
                << "] - [ "
                << static_cast<cellPtrT>(_rootPtr)->xCenter + 0.5*rootSize
                << ","
                << static_cast<cellPtrT>(_rootPtr)->yCenter + 0.5*rootSize
                << ","
                << static_cast<cellPtrT>(_rootPtr)->zCenter + 0.5*rootSize
                << "] )";
  Logger.flushStream();
  Logger.destroy();
};

PartsTooClose::~PartsTooClose()
{
};


};

#endif
