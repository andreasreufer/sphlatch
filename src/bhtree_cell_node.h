#ifndef BHTREE_CELL_NODE_H
#define BHTREE_CELL_NODE_H

/*
 *  bhtree_cell_node.h
 *  
 *
 *  Created by Andreas Reufer on 08.02.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */


#include "typedefs.h"
#include "bhtree_generic_node.h"

namespace sphlatch {
/// 
/// \brief generic cell node
///
struct genericCellNode : public genericNode {
  typedef genericNode* genericNodePtrType;

///
/// pointers to children
///
  genericNodePtrType child[8];

///
/// center and size of the cell
///
  fType xCenter, yCenter, zCenter;
  fType cellSize;
};

};

#endif
