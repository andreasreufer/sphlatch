#ifndef BHTREE_PARTICLE_PROXY_H
#define BHTREE_PARTICLE_PROXY_H

/*
 *  bhtree_particle_proxy.h
 *
 *
 *  Created by Andreas Reufer on 15.11.07.
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"

#include "particle_node.h"

namespace sphlatch {
struct particleProxy {

  typedef particleProxy*      particleProxyPtrType;
  typedef particleNode*       particleNodePtrType;

  matrixPtrType matrixPtr;
  size_t rowIndex;

  particleProxy(particleProxymatrixPtrType _matrixPtr, size_t const _rowIndex)
  {
    matrixPtr = _matrixPtr;
    rowIndex = _rowIndex;
  }

  particleProxy(void)
  {
  }

  particleProxy* operator*()
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

