#ifndef OCTREE_H
#define OCTREE_H

/*
 *  octree.h
 *  
 *
 *  Created by Andreas Reufer on 12.11.07.
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>

#include "typedefs.h"

#include "nodeproxy.h"
typedef NodeProxy	NodeProxyType;
typedef NodeProxy*	NodeProxyPtrType;
typedef NodeProxy&	NodeProxyRefType;

#include "octreenode.h"
typedef GenericOctNode<NodeProxyPtrType>* GenericOctNodePtr;

enum ParticleIndex { PID, X, Y, Z, M, AX, AY, AZ, PSIZE };
enum MonopolIndex  { CID, CX, CY, CZ, Q000, MSIZE};

//template <typename T>
class OctTree {
	
	public:
		/** \brief constructor:
		* - allocate and setup RootNode and point cursor to it
		* - set counters
		* - instantiate recursors
		* - build toptree
		* - set up CellMultipoles matrix
		*/
		OctTree(void) {
			//RootPtr = new GenericOctNode<T>;
			RootPtr = new GenericOctNode<NodeProxyPtrType>;
			
			RootPtr->parent = NULL;
			RootPtr->child[0] = NULL;
			RootPtr->child[1] = NULL;
			RootPtr->child[2] = NULL;
			RootPtr->child[3] = NULL;
			RootPtr->child[4] = NULL;
			RootPtr->child[5] = NULL;
			RootPtr->child[6] = NULL;
			RootPtr->child[7] = NULL;
		
			RootPtr->payload = NULL;
			
			// get this from costzone!
			RootPtr->xCenter = 0.5;
			RootPtr->yCenter = 0.5;
			RootPtr->zCenter = 0.5;
			RootPtr->cellSize = 1.;
			
			RootPtr->isParticle		= false;
			RootPtr->isEmpty		= true;
			RootPtr->isGravitating	= false;
			RootPtr->depth = 0;
			
			CurNodePtr = RootPtr;
			
			nodeCounter	= 1;
			particleCounter	= 0;
			
			toptreeDepth = 0; // get this from costzone and MAC!
			//thetaMAC = 1.;
			thetaMAC = 0.4;
			
			buildToptreeRecursor();
		};
		
		/**
		* destructor:
		*  - postorder recurse node deletion
		*/
		~OctTree(void) {
			goRoot();
			empty();
			goRoot();
			delete CurNodePtr; // Seppuku!
		}
		
	private:
		GenericOctNodePtr CurNodePtr, RootPtr;
			
		/**
		* variables
		*/
		size_t nodeCounter, particleCounter, toptreeDepth;
		size_t debugCounter;	// just to fool around
		
		matrixType				cellData;
		std::vector<NodeProxy>	cellProxies;
		
		   
		// little private helpers
	private:
		/**
		* go up one level
		*/
		inline void goUp() {
			CurNodePtr = CurNodePtr->parent;
		}

		/**
		* go to child
		*/
		inline void goChild(const size_t _n) {
			CurNodePtr = CurNodePtr->child[_n];
		}

		/**
		* go to root
		*/
		inline void goRoot() {
			CurNodePtr = RootPtr;
		}

		/**
		* setup coordinates and cell size as child _n
		*/
		void setupCoordinates(const size_t _n) {
			if ( CurNodePtr->parent != NULL ) {
				CurNodePtr->cellSize = CurNodePtr->parent->cellSize / 2.;
				CurNodePtr->xCenter = ( ( _n      )% 2 ) ?
					CurNodePtr->parent->xCenter + 0.5*CurNodePtr->cellSize :
					CurNodePtr->parent->xCenter - 0.5*CurNodePtr->cellSize ; 
				CurNodePtr->yCenter = ( ( _n >> 1 )% 2 ) ?
					CurNodePtr->parent->yCenter + 0.5*CurNodePtr->cellSize :
					CurNodePtr->parent->yCenter - 0.5*CurNodePtr->cellSize ; 
				CurNodePtr->zCenter = ( ( _n >> 2 )% 2 ) ?
					CurNodePtr->parent->zCenter + 0.5*CurNodePtr->cellSize :
					CurNodePtr->parent->zCenter - 0.5*CurNodePtr->cellSize ; 
			}
		}

		/**
		* create a new child in octant _n
		*/
		void newChild(const size_t _n) {
			if ( CurNodePtr->child[_n] == NULL ) {
				//GenericOctNodePtr NewNodePtr = new GenericOctNode<T>;
				GenericOctNodePtr NewNodePtr = 
					new GenericOctNode<NodeProxyPtrType>;

				NewNodePtr->parent = CurNodePtr;
				CurNodePtr->child[_n] = NewNodePtr;
				
				goChild(_n);

				CurNodePtr->child[0] = NULL;
				CurNodePtr->child[1] = NULL;
				CurNodePtr->child[2] = NULL;
				CurNodePtr->child[3] = NULL;
				CurNodePtr->child[4] = NULL;
				CurNodePtr->child[5] = NULL;
				CurNodePtr->child[6] = NULL;
				CurNodePtr->child[7] = NULL;
				
				CurNodePtr->depth = CurNodePtr->parent->depth + 1;
				
				setupCoordinates(_n);
				
				goUp();
			}
		}
		
		/**
		* convert current node to cell, if current node is a particle
		*/
		void converttoCell(const size_t _n) {
			if ( CurNodePtr->isParticle ) {
				setupCoordinates(_n);
				
				CurNodePtr->payload = NULL;
				
				CurNodePtr->isParticle = false;
				CurNodePtr->isEmpty    = true;
				
				nodeCounter++;
			}
		};
		// end of little helpers
		
		
		// top tree stuff
	private:
		/**
		* recursor to build toptree
		*/
		void buildToptreeRecursor(void) {
			if ( belowToptreeStop() ) {} else {
				makeEmptyCells();
				for (size_t i = 0; i < 8; i++) { // try without loop
					if ( CurNodePtr->child[i] != NULL ) {
						goChild(i);
						buildToptreeRecursor();
						goUp();
					}
				}
			}
		};	
		/**
		* make empty cell nodes for the toptree
		*/
		void makeEmptyCells(void) {
			for (size_t i = 0; i < 8; i++) {
				newChild(i);
				goChild(i);
				CurNodePtr->payload = NULL;
				CurNodePtr->isParticle		= false;
				CurNodePtr->isEmpty			= true;
				CurNodePtr->isGravitating	= false;
				nodeCounter++;
				goUp();
			}
		}
		/**
		* stop recursion if:
		* - depth of current node is below toptreeDepth
		*/
		bool belowToptreeStop(void) {
			return ( CurNodePtr->depth >= toptreeDepth );
		};
		// end of top tree stuff
	
		
		// insertParticle() stuff
	public:
		/**
		* method to insert particle:
		*  - go to root
		*  - call the insertion recursor
		*/		
		void insertParticle(NodeProxyType* newPayload, 
				nodetypeType newType) {
			goRoot();
			insertParticleRecursor(newPayload, newType);
			particleCounter++;
		}
	private:		
		/**
		* recursor for inserting a new particle:
		* try to insert as child of current
		* node. if child is
		*	- empty, insert particle. we're done.
		*	- a node, go to child and call recursor
		*	- a particle, disconnect particle and
		*	  call recursor for existing and new
		*	  particle.
		*/
		void insertParticleRecursor(NodeProxyType* _newPayload,
				bool _newIsGravitating) {
			size_t targetOctant = 0;
			targetOctant += ( (*_newPayload)(X) < CurNodePtr->xCenter ) ? 0 : 1;
			targetOctant += ( (*_newPayload)(Y) < CurNodePtr->yCenter ) ? 0 : 2;
			targetOctant += ( (*_newPayload)(Z) < CurNodePtr->zCenter ) ? 0 : 4;

			/**
			* If targeted child is empty, place the particle there
			**/
			if ( CurNodePtr->child[targetOctant] == NULL ) {
				CurNodePtr->isEmpty			= false;
				
				newChild(targetOctant);
				goChild(targetOctant);

				CurNodePtr->payload			= _newPayload;
				CurNodePtr->isParticle		= true;
				//CurNodePtr->isEmpty			= false;
				CurNodePtr->isEmpty			= true;
				CurNodePtr->isGravitating	= _newIsGravitating;
				
				/* particle saves its position to node directly */
				CurNodePtr->xCom = (*_newPayload)(X);
				CurNodePtr->yCom = (*_newPayload)(Y);
				CurNodePtr->zCom = (*_newPayload)(Z);
				CurNodePtr->q000 = (*_newPayload)(M);				
				
				/**
				* don't forget to wire the nodePtr of the 
				* NodeProxy back to the node
				*/
				CurNodePtr->payload->nodePtr = CurNodePtr;
				
				goUp();
			}
		
			/**
			* ... or if existing child is a node, then try to place the particle
			* as a child of this node
			*/
			else if ( ! CurNodePtr->child[targetOctant]->isParticle ) {
				CurNodePtr->isEmpty = false;
				goChild(targetOctant);
				insertParticleRecursor(_newPayload, _newIsGravitating);
				goUp();
			}
			
			/**
			* ... or if existing child is a particle (ghost/nonghost), then
			* replace it by a new node and try to place the existing two
			* particles as childs of this node
			*/
			else if ( CurNodePtr->child[targetOctant]->isParticle ) {
				/**
				* goto child, save resident particle
				* and convert it to a node
				*/
				goChild(targetOctant);
				
				NodeProxyType*	residentPayload = CurNodePtr->payload;
				bool	residentGravitating	= CurNodePtr->isGravitating;
				
				converttoCell(targetOctant);

				insertParticleRecursor(residentPayload, residentGravitating);
				insertParticleRecursor(_newPayload, _newIsGravitating);
			} else {
				/* this never happens */
			}
		}		
		// end of insertParticle() stuff
		
		
		// calcMultipoles() stuff
	public:
		/**
		* calculate multipoles:
		*  - prepare cellData matrix and cellProxies vector
		*  - connect every node to a matrix row
		*  - go to root
		*  - call the multipole recursor
		*  - exchange toptrees
		*/
		void calcMultipoles(void) {
			cellData.resize(nodeCounter, MSIZE);
			cellProxies.resize(nodeCounter);

			for (size_t i = 0; i < nodeCounter; i++) {
				( cellProxies[i] ).setup(&cellData, i);
				cellData(i, CID) = i;
			}
			
			nodeCounter = 0;
			goRoot();
			connectCellDataRecursor();
			
			debugCounter = 0;
			goRoot();
			calcMultipoleRecursor();
#ifdef OOSPH_MPI
			//(Exchange)
#endif
		}
	
	private:
		/**
		* recursor for multipole calculation
		*/
		void calcMultipoleRecursor(void) {
			if ( calcMultipoleStop() ) {} else {
				for (size_t i = 0; i < 8; i++) { // try without loop
					if ( CurNodePtr->child[i] != NULL ) {
						goChild(i);
						calcMultipoleRecursor();
						goUp();
					}
				}
				calcMultipole();
			}
		};	
		
		/**
		* stop recursion if:
		* - current node is empty
		* - current node is a particle
		*/
		bool calcMultipoleStop(void) {
			return ( CurNodePtr->isEmpty || CurNodePtr->isParticle );
		};
		
		/**
		* recursive function to connect
		* cells to a matrix row
		*/
		void connectCellDataRecursor(void) {
			connectCellData();
			for (size_t i = 0; i < 8; i++) { // try without loop
				if ( CurNodePtr->child[i] != NULL ) {
					goChild(i);
					connectCellDataRecursor();
					goUp();
				}
			}
		}

		/**
		* connect cell to a matrix row
		*/
		void connectCellData() {
			if ( ! (CurNodePtr->isParticle) ) {
				CurNodePtr->payload = *(cellProxies[nodeCounter]);
				nodeCounter++;
			}
		}

		/**
		* calculate multipole from children
		* todo: include isGravitating flag
		*/
		valueType monopolCM, monopolCXM, monopolCYM, monopolCZM;
		void calcMultipole() {
			/*if (! CurNodePtr->isParticle ) {*/	// save this check by making
				if (! CurNodePtr->isEmpty ) {	//  particles empty
					monopolCM = 0.;
					monopolCXM = 0.;
					monopolCYM = 0.;
					monopolCZM = 0.;
					for (size_t i = 0; i < 8; i++) {
						if ( CurNodePtr->child[i] != NULL &&
							CurNodePtr->isGravitating ) {
							goChild(i);
							monopolCM += CurNodePtr->q000;
							monopolCXM += (CurNodePtr->xCom)*
											(CurNodePtr->q000);
							monopolCYM += (CurNodePtr->yCom)*
											(CurNodePtr->q000);
							monopolCZM += (CurNodePtr->zCom)*
											(CurNodePtr->q000);
							goUp();
						}
					}
										
					/* copy data to node itself */
					CurNodePtr->q000 = monopolCM;
					CurNodePtr->xCom = monopolCXM / monopolCM;
					CurNodePtr->yCom = monopolCYM / monopolCM;
					CurNodePtr->zCom = monopolCZM / monopolCM;					
				}
			//}
		}			
		// end of multipole stuff

		
		// calcGravity() stuff
	private:
		NodeProxyType* curGravNodeProxy;
		valueType curGravParticleX, curGravParticleY, curGravParticleZ;
		valueType curGravParticleAX, curGravParticleAY, curGravParticleAZ;
		valueType thetaMAC;
		size_t calcGravityCellsCounter, calcGravityPartsCounter;
	public:
		/**
		* calculate gravitation for a particle:
		*  - load current particle data
		*  - go to root
		*  - call the gravity calculation recursor
		*  - write back resulting acceleration
		*/
		void calcGravity(NodeProxyType* _curParticle) {
			curGravNodeProxy = _curParticle;
						
			curGravParticleX = (*curGravNodeProxy)(X);
			curGravParticleY = (*curGravNodeProxy)(Y);
			curGravParticleZ = (*curGravNodeProxy)(Z);
			curGravParticleAX = 0.;
			curGravParticleAY = 0.;
			curGravParticleAZ = 0.;
#ifdef TREEPROFILE
			calcGravityPartsCounter = 0;
			calcGravityCellsCounter = 0;
#endif
			
			/**
			* trick: hide the current particle by letting it look like
			* an empty cell node, so that it doesn't gravitate with
			* itself.
			*/
			curGravNodeProxy->nodePtr->isParticle	= false;
			curGravNodeProxy->nodePtr->isEmpty		= true;
			
			calcGravityPartsCounter = 0;
			calcGravityCellsCounter = 0;
			
			goRoot();
			calcGravityRecursor();
			
			(*curGravNodeProxy)(AX) += curGravParticleAX;
			(*curGravNodeProxy)(AY) += curGravParticleAY;
			(*curGravNodeProxy)(AZ) += curGravParticleAZ;
			
			/**
			* de-trick
			*/
			curGravNodeProxy->nodePtr->isParticle	= true;
			curGravNodeProxy->nodePtr->isEmpty		= false;			
		}
		
	private:
		void calcGravityRecursor(void) {
			if ( CurNodePtr->isParticle ) {
				calcGravParticle();
			} else
			if ( CurNodePtr->isEmpty ) {
			} else
			if ( calcGravMAC() ) {
				calcGravCell();					
			} else {
				for (size_t i = 0; i < 8; i++) { // try without loop
					if ( CurNodePtr->child[i] != NULL ) {
						goChild(i);
						calcGravityRecursor();
						goUp();
					}
				}
			}
		}
		/**
		* stop recursion if:
		* - current node is empty
		* - MAC is fulfilled
		*/
		valueType cellPartDist;
		bool calcGravMAC(void) {
			cellPartDist = sqrt( (CurNodePtr->xCenter - curGravParticleX)*
				(CurNodePtr->xCenter - curGravParticleX) +
				(CurNodePtr->yCenter - curGravParticleY)*
				(CurNodePtr->yCenter - curGravParticleY) + 
				(CurNodePtr->zCenter - curGravParticleZ)*
				(CurNodePtr->zCenter - curGravParticleZ)
				);
			return ( ( ( CurNodePtr->cellSize) / cellPartDist ) < thetaMAC );
		}

		/**
		* calculate acceleration due to a particle
		* todo: add grav-const
		*/
		valueType partGravPartnerX, partGravPartnerY, partGravPartnerZ, 
				  partGravPartnerM;
		void calcGravParticle() {
#ifdef TREEPROFILE
			calcGravityPartsCounter++;
#endif
			partGravPartnerX = CurNodePtr->xCom;
			partGravPartnerY = CurNodePtr->yCom;
			partGravPartnerZ = CurNodePtr->zCom;
			partGravPartnerM = CurNodePtr->q000;
			
			cellPartDist = sqrt(	(partGravPartnerX - curGravParticleX)*
				(partGravPartnerX - curGravParticleX) +
				(partGravPartnerY - curGravParticleY)*
				(partGravPartnerY - curGravParticleY) + 
				(partGravPartnerZ - curGravParticleZ)*
				(partGravPartnerZ - curGravParticleZ)
				); // include softening here
									
			cellPartDistPow3 = cellPartDist*cellPartDist*cellPartDist;
			
			curGravParticleAX -=  partGravPartnerM *
				( curGravParticleX - partGravPartnerX ) / 
				cellPartDistPow3;
			curGravParticleAY -=  partGravPartnerM *
				( curGravParticleY - partGravPartnerY ) / 
				cellPartDistPow3;
			curGravParticleAZ -=  partGravPartnerM *
				( curGravParticleZ - partGravPartnerZ ) / 
				cellPartDistPow3;
		}
		
		/**
		* calculate acceleration due to a cell
		* todo: add grav-const
		*/
		valueType cellPartDistPow3;
		void calcGravCell() {
#ifdef TREEPROFILE
			calcGravityCellsCounter++;
#endif
			curGravParticleAX -= ( CurNodePtr->q000 ) * 
				( curGravParticleX - CurNodePtr->xCom ) /
				cellPartDistPow3;
			curGravParticleAY -= ( CurNodePtr->q000 ) * 
				( curGravParticleY - CurNodePtr->yCom ) /
				cellPartDistPow3;
			curGravParticleAX -= ( CurNodePtr->q000 ) * 
				( curGravParticleZ - CurNodePtr->zCom ) /
				cellPartDistPow3;
			
		}
		// end of calcGravity() stuff
		
		// neighbour search stuff
	public:
		/**
		* find neighbours:
		*  - load current particle data
		*  - go to current particle node
		*  - go up until every neighbour is within current cell
		*  - call the neighbour search recursor
		*  - brute force sort out non-neighbours
		*  - return neighbours
		*/
		/*somecontainer_type findNeighbours(
				const GenericOctNodePtr _curParticle, 
				const valueRefType _search_radius) {
			CurNodePtr = _curParticle;
			size_t topDepth = <crazyformula>
			while ( CurNodePtr->depth > topDepth ) {
				goUp();
			}
			neighbourFindRecursor(_curParticle);
			return neighbours;
		}*/
	private:
		// end of neighbour search stuff


		// empty() stuff
		/**
		* deletes the current subtree
		*/
	private:	
		void empty(void) {
			emptyRecursor();
		}
		
		/**
		* recursor for emptying the
		* tree
		*/
		void emptyRecursor() {
			for (size_t i = 0; i < 8; i++) { // try without loop
				if ( CurNodePtr->child[i] != NULL ) {
					goChild(i);
					emptyRecursor();
					goUp();
					
					delete CurNodePtr->child[i];
					CurNodePtr->child[i] = NULL;
				}
			}
		};
		// end of empty() stuff
};

#endif
