#ifndef PARTICLEPROXY_H
#define PARTICLEPROXY_H

/*
 *  particleproxy.h
 *  
 *
 *  Created by Andreas Reufer on 15.11.07.
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"

#include "octreenode.h"

struct ParticleProxy {
	typedef ParticleProxy*	ParticleProxyTypePtr;
	GenericOctNode<ParticleProxyTypePtr>*	nodePtr;
	
	matrixPtrType	matrixPtr;
	size_t			rowIndex;

	ParticleProxy(matrixPtrType _matrixPtr, size_t const _rowIndex) {
		matrixPtr = _matrixPtr;
		rowIndex = _rowIndex;
	}

	ParticleProxy(void) {}

	ParticleProxy* operator*() {
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