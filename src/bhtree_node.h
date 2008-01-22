#ifndef BHTREE_NODE_H
#define BHTREE_NODE_H

/*
 *  bhtree_node.h
 *  
 *
 *  Created by Andreas Reufer on 15.11.07.
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"

namespace sphlatch {

/** \brief tree node struct with an
*           arbitrary payload
*/

template <typename T>
struct GenericOctNode {
	typedef GenericOctNode* GenericOctNodePtr;
	/**
	* pointers to parent and children
	*/
	GenericOctNodePtr parent;
	GenericOctNodePtr child[8];
	
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
	* in case of a cell node, the coordinates are
	* set to the center of the cell.  in case of 
	* particle the coordinates of the particle are
	* saved for caching reasons.
	*/
	valueType xCenter;
	valueType yCenter;
	valueType zCenter;
	valueType cellSize;

    /**
	* payload of the node
	*/
	T payload;
	
	valueType xCom, yCom, zCom,	q000;
	valueType ax, ay, az;

	/**
	* isParticle:    is a particle
	*/
	bool isParticle		:1;	
	 
	/**
	* isEmpty:		 a cell node whose subtrees do not contain particles
	*                on any costzone domain. 
	*/
	bool isEmpty		:1;
	
	/**
	* isLocal:       determines whether a particle is local or non-
    *                local (ghost). 
	*/
	bool isLocal	:1;
	
	/**
	* pointer operator
	*/
	GenericOctNode*  operator*() {
		return this;
	}
};

};

#endif

