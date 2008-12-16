#ifndef MONOPOLE_NODE_H
#define MONOPOLE_NODE_H

/*
 *  monopole_node.h
 *
 *
 *  Created by Andreas Reufer on 06.02.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"
#include "bhtree_cell_node.h"

namespace sphlatch {
///
/// \brief monopole cell node
///
struct monopoleCellNode : public genericCellNode {
///
/// center of mass and mass
///
  fType xCom, yCom, zCom;
  fType mass;
};
};

#endif

