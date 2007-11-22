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

#include "particleproxy.h"
typedef ParticleProxy	ParticleProxyType;
typedef ParticleProxy*	ParticleProxyPtrType;
typedef ParticleProxy&	ParticleProxyRefType;

#include "octreenode.h"
typedef GenericOctNode<ParticleProxyPtrType>* GenericOctNodePtr;

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
			RootPtr = new GenericOctNode<ParticleProxyPtrType>;
			
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
			
			toptreeDepth = 1; // get this from costzone and MAC!
			//thetaMAC = 1.;
			thetaMAC = 0.7;
			
			buildToptreeRecursor();
		};
		
		/**
		* destructor:
		*  - postorder recurse node deletion
		*/
		~OctTree(void) {}
	
		/**
		* method to insert particle:
		*  - go to root
		*  - call the insertion recursor
		*/		
		void insertParticle(ParticleProxyType* newPayload, 
				nodetypeType newType) {
			goRoot();
			insertParticleRecursor(newPayload, newType);
			particleCounter++;
		}

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
			//std::cout << nodeCounter << " nodes \n";
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
			/*std::cout << "empty: " << debugCounter << "\n";
			
			std::cout << (*(RootPtr->payload))(CX) << "   "
				<< (*(RootPtr->payload))(CY) << "   "
				<< (*(RootPtr->payload))(CZ) << "   "
				<< (*(RootPtr->payload))(Q000) << "\n";*/
			
#ifdef OOSPH_MPI
			//(Exchange)
#endif
		}
		 
		/**
		* calculate gravitation for a particle:
		*  - load current particle data
		*  - go to root
		*  - call the gravity calculation recursor
		*  - write back resulting acceleration
		*/
		ParticleProxyType* curGravParticleProxy;
		valueType curGravParticleX, curGravParticleY, curGravParticleZ;
		valueType curGravParticleAX, curGravParticleAY, curGravParticleAZ;
		valueType thetaMAC;
		
		size_t calcGravityCellsCounter, calcGravityPartsCounter;
		void calcGravity(ParticleProxyType* _curParticle) {
			curGravParticleProxy = _curParticle;
						
			curGravParticleX = (*curGravParticleProxy)(X);
			curGravParticleY = (*curGravParticleProxy)(Y);
			curGravParticleZ = (*curGravParticleProxy)(Z);
			curGravParticleAX = 0.;
			curGravParticleAY = 0.;
			curGravParticleAZ = 0.;
						
			// do something so that _curParticle
			// is excluded from the recursion
			
			/**
			* trick: hide the current particle by letting it look like
			* an empty cell node, so that it doesn't gravitate with
			* itself.
			*/
			curGravParticleProxy->nodePtr->isParticle	= false;
			curGravParticleProxy->nodePtr->isEmpty		= true;
			
			calcGravityPartsCounter = 0;
			calcGravityCellsCounter = 0;
			
			goRoot();
			calcGravityRecursor();
			
			/*std::cout << "\nstats: " << calcGravityPartsCounter
				<< " particles and "
				<< calcGravityCellsCounter << " cells\n";*/
			
			(*curGravParticleProxy)(AX) += curGravParticleAX;
			(*curGravParticleProxy)(AY) += curGravParticleAY;
			(*curGravParticleProxy)(AZ) += curGravParticleAZ;
			
			/**
			* de-trick
			*/
			curGravParticleProxy->nodePtr->isParticle	= true;
			curGravParticleProxy->nodePtr->isEmpty		= false;			
		}
	
		  
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

		/**
		* stop recursion if:
		* - when we've arrived in Neverland :-)
		*/
		bool neverStop(void) {
			return false;
		};
		
		/**
		* stop recursion if:
		* - depth of current node is below toptreeDepth
		*/
		bool belowToptreeStop(void) {
			return ( CurNodePtr->depth >= toptreeDepth );
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
		* stop recursion if:
		* - current node is empty
		* - MAC is fulfilled
		*/
		/*valueTye theta;
		bool calcGravStop(void) {
			if ( CurNodePtr->isEmpty || CurNodePtr->isParticle ) {
				return true;
			} else
			if ( 
				
				||  CurNodePtr->isParticle
			return ( CurNodePtr->isEmpty || CurNodePtr->isParticle );
			//		insert MAC
		}*/
		
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
		* do nothing :-)
		*/
		void doNothing(void) {}
		
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
			if (! CurNodePtr->isParticle ) {	// save this check by making
				if (! CurNodePtr->isEmpty ) {	//  particles empty
					monopolCM = 0.;
					monopolCXM = 0.;
					monopolCYM = 0.;
					monopolCZM = 0.;
					for (size_t i = 0; i < 8; i++) {
						if ( CurNodePtr->child[i] != NULL ) {
							goChild(i);
							if ( CurNodePtr->isParticle ) {
								monopolCM  += (*(CurNodePtr->payload))(M);
								monopolCXM += (*(CurNodePtr->payload))(X) *
									(*(CurNodePtr->payload))(M);
								monopolCYM += (*(CurNodePtr->payload))(Y) *
									(*(CurNodePtr->payload))(M);
								monopolCZM += (*(CurNodePtr->payload))(Z) *
									(*(CurNodePtr->payload))(M);
								} else
							if (! CurNodePtr->isEmpty ) {
								monopolCM  += (*(CurNodePtr->payload))(Q000);
								monopolCXM += (*(CurNodePtr->payload))(CX) *
									(*(CurNodePtr->payload))(Q000);
								monopolCYM += (*(CurNodePtr->payload))(CY) *
									(*(CurNodePtr->payload))(Q000);
								monopolCZM += (*(CurNodePtr->payload))(CZ) *
									(*(CurNodePtr->payload))(Q000);
							} else {}
							goUp();
						}
					}
					(*(CurNodePtr->payload))(Q000) = monopolCM;
					(*(CurNodePtr->payload))(CX) = monopolCXM / monopolCM;
					(*(CurNodePtr->payload))(CY) = monopolCYM / monopolCM;
					(*(CurNodePtr->payload))(CZ) = monopolCZM / monopolCM;
				}
			}
		}			
		
		/**
		* calculate acceleration due to a particle
		* todo: add grav-const
		*/
		
		valueType partGravPartnerX, partGravPartnerY, partGravPartnerZ, 
				  partGravPartnerM;
		void calcGravParticle() {
			//std::cout << CurNodePtr->depth << "p ";
			//calcGravityPartsCounter++;
			partGravPartnerX = (*(CurNodePtr->payload))(X);
			partGravPartnerY = (*(CurNodePtr->payload))(Y);
			partGravPartnerZ = (*(CurNodePtr->payload))(Z);
			partGravPartnerM = (*(CurNodePtr->payload))(M);
			
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
			//std::cout << CurNodePtr->depth << "c ";
			//calcGravityCellsCounter++;
			cellPartDistPow3 = cellPartDist*cellPartDist*cellPartDist;
			curGravParticleAX -=  (*(CurNodePtr->payload))(Q000) *
				( curGravParticleX - (*(CurNodePtr->payload))(CX) ) / 
				cellPartDistPow3;
			curGravParticleAY -=  (*(CurNodePtr->payload))(Q000) *
				( curGravParticleY - (*(CurNodePtr->payload))(CY) ) / 
				cellPartDistPow3;
			curGravParticleAZ -=  (*(CurNodePtr->payload))(Q000) *
				( curGravParticleZ - (*(CurNodePtr->payload))(CZ) ) / 
				cellPartDistPow3;
		}
		
		/**
		* go up one level if there is a parent
		*/
		inline void goUp() {
			CurNodePtr = CurNodePtr->parent;
		}

		/**
		* go to child if it exists
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
					new GenericOctNode<ParticleProxyPtrType>;

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
		* recursor for inserting a new particle:
		* try to insert as child of current
		* node. if child is
		*	- empty, insert particle. we're done.
		*	- a node, go to child and call recursor
		*	- a particle, disconnect particle and
		*	  call recursor for existing and new
		*	  particle.
		*/
		void insertParticleRecursor(ParticleProxyType* _newPayload,
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
				CurNodePtr->isEmpty			= false;
				CurNodePtr->isGravitating	= _newIsGravitating;
				
				/**
				* don't forget to wire the nodePtr of the 
				* ParticleProxy back to the node
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
				
				ParticleProxyType*	residentPayload = CurNodePtr->payload;
				bool	residentGravitating	= CurNodePtr->isGravitating;
				
				converttoCell(targetOctant);

				insertParticleRecursor(residentPayload, residentGravitating);
				insertParticleRecursor(_newPayload, _newIsGravitating);
			} else {
				/* this never happens */
			}
		}
		
		/**
		* */
		GenericOctNodePtr CurNodePtr, RootPtr;
			
		/**
		* variables
		*/
		size_t nodeCounter, particleCounter, toptreeDepth;
		size_t debugCounter;	// just to fool around
		
		matrixType					cellData;
		std::vector<ParticleProxy>	cellProxies;
};

#endif
