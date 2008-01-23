#ifndef BHTREE_NODE_PROXY_H
#define BHTREE_NODE_PROXY_H

/*
 *  bhtree_node_proxy.h
 *
 *
 *  Created by Andreas Reufer on 15.11.07.
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"

#include "bhtree_node.h"

namespace sphlatch {
struct NodeProxy {
  typedef NodeProxy*      NodeProxyTypePtr;
  GenericOctNode<NodeProxyTypePtr>*       nodePtr;

  matrixPtrType matrixPtr;
  size_t rowIndex;

  NodeProxy(matrixPtrType _matrixPtr, size_t const _rowIndex)
  {
    matrixPtr = _matrixPtr;
    rowIndex = _rowIndex;
  }

  NodeProxy(void)
  {
  }

  NodeProxy* operator*()
  {
    return this;
  }

  valueRefType operator ()(const size_t &j)
  {
    return (*matrixPtr)(rowIndex, j);
  }

  void setup(matrixPtrType _matrixPtr, size_t const _rowIndex)
  {
    matrixPtr = _matrixPtr;
    rowIndex = _rowIndex;
  }
};
};

#endif

