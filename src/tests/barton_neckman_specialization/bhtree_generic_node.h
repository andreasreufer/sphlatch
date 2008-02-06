#ifndef BHTREE_GENERIC_NODE_H
#define BHTREE_GENERIC_NODE_H

/*
 *  bhtree_generic_node.h
 *
 *
 *  Created by Andreas Reufer on 06.02.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"

namespace sphlatch {
/** \brief tree node struct with an
 *           arbitrary payload
 */

struct genericNode {
  typedef genericNode* genericNodePtrType;

  /**
   * pointer to parent
   */
  genericNodePtrType parent;

  /**
   * the root node has depth = 0, its
   * children have depth = 1 and so on
   */
  size_t depth;

  /**
   * identifier of the node which is casted from
   * the ID row for particles. should be unique in
   * each local tree.
   */
  identType ident;

  /**
   * payload of the node
   */
  T payload;

  valueType xCom, yCom, zCom, q000;

  /**
   * isParticle:    is a particle
   */
  bool isParticle : 1;

  /**
   * isEmpty:		 a cell node whose subtrees do not contain particles
   *                on any costzone domain.
   */
  bool isEmpty    : 1;

  /**
   * isLocal:       determines whether a particle is local or non-
   *                local (ghost).
   */
  bool isLocal    : 1;

  /**
   * pointer operator
   */
  genericNode*  operator*() {
    return this;
  }
};

struct genericCellNode : public genericNode {
  typedef genericNode* genericNodePtrType;

  /**
   * pointers to children
   */
  genericNodePtrType child[8];

  /**
   * center and size of the cell
   */
  valueType xCenter, yCenter, zCenter;
  valueType cellSize;
};
};

#endif

