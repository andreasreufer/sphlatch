#ifndef QUADRUPOLE_NODE_H
#define QUADRUPOLE_NODE_H

/*
 *  quadrupole_node.h
 *
 *
 *  Created by Andreas Reufer on 06.02.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"
#include "monopole_node.h"

namespace sphlatch {
/** \brief quadrupole cell node
 */

struct quadrupoleCellNode : public monopoleCellNode {
  /**
   * quadrupole terms
   */
  valueType q11, q22, q33, q12, q13, q23;
};
};

#endif

