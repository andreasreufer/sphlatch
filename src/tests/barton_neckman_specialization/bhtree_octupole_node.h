#ifndef OCTUPOLE_NODE_H
#define OCTUPOLE_NODE_H

/*
 *  octupole_node.h
 *
 *
 *  Created by Andreas Reufer on 06.02.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"
#include "bhtree_quadrupole_node.h"

namespace sphlatch {
///
/// \brief octupole cell node
///
struct octupoleCellNode : public quadrupoleCellNode {
///
/// octupole moments
///
  valueType s11, s22, s33, s12, s21, s13, s31, s23, s32, s123;
};
};

#endif

