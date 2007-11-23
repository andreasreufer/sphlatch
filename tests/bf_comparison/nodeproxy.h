#ifndef NODEPROXY_H
#define NODEPROXY_H

/*
 *  nodeproxy.h
 *  
 *
 *  Created by Andreas Reufer on 15.11.07.
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"

#include "octreenode.h"
struct NodeProxy {
	typedef NodeProxy*	NodeProxyTypePtr;
	GenericOctNode<NodeProxyTypePtr>*	nodePtr;
	
	matrixPtrType	matrixPtr;
	size_t			rowIndex;

	NodeProxy(matrixPtrType _matrixPtr, size_t const _rowIndex) {
		matrixPtr = _matrixPtr;
		rowIndex = _rowIndex;
	}

	NodeProxy(void) {}

	NodeProxy* operator*() {
		return this;
	}
	
	valueRefType operator ()(const size_t &j) {
		return  (*matrixPtr)(rowIndex, j);
	}

	void setup(matrixPtrType _matrixPtr, size_t const _rowIndex) {
		matrixPtr = _matrixPtr;
		rowIndex = _rowIndex;
	}
};

#endif