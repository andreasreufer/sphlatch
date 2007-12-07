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
#include <stack>
#include <queue>

#ifdef OOSPH_MPI
#include <mpi.h>
#endif

#include <fstream>

#include "typedefs.h"

#include "nodeproxy.h"
typedef NodeProxy	NodeProxyType;
typedef NodeProxy*	NodeProxyPtrType;
typedef NodeProxy&	NodeProxyRefType;

#include "octreenode.h"
typedef GenericOctNode<NodeProxyPtrType>* GenericOctNodePtr;

enum ParticleIndex { PID, X, Y, Z, M, AX, AY, AZ, PSIZE };
enum MonopolIndex  { CID, CX, CY, CZ, Q000, MSIZE};
//enum QuadrupolIndex  { CID, CX, CY, CZ, Q000, Q001, Q002, Q010, Q011, MSIZE};

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
			RootPtr->isLocal    	= true;
			RootPtr->depth = 0;
			
			CurNodePtr = RootPtr;
			
			cellCounter	= 1;  // we now have the root cell and no particles
			partCounter	= 0;
			
			toptreeDepth = 6; // get this from costzone and MAC!
			thetaMAC = 1.;
			//thetaMAC = 0.4;
			
			buildToptreeRecursor();
            noToptreeCells = cellCounter;
            cellIsFilled.resize(noToptreeCells);
		};
		
		/**
		* destructor:
		*  - postorder recurse node deletion
		*/
		~OctTree(void) {
			goRoot();
			empty();
			delete CurNodePtr; // Seppuku!
		}
		
	private:
		GenericOctNodePtr CurNodePtr, RootPtr;
			
		/**
		* variables
		*/
		size_t cellCounter, partCounter, toptreeDepth, noToptreeCells;
		size_t debugCounter;	// just to fool around
		
		matrixType				cellData;
		std::vector<NodeProxy>	cellProxies;
        bitsetType cellIsFilled; // bool vector for toptree
		
		   
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
				
				cellCounter++;
			}
		};
		// end of little helpers
		
		
		// top tree building stuff
	private:
		/**
		* recursor to build toptree
		*/
		void buildToptreeRecursor(void) {
			if ( atToptreeStop() ) {} else {
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
				CurNodePtr->isLocal	= false;
				cellCounter++;
				goUp();
			}
		}
		/**
		* stop recursion if:
		* - depth of current node is below toptreeDepth
		*/
		bool atToptreeStop(void) {
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
			partCounter++;
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
				bool _newIsLocal) {
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
				CurNodePtr->isLocal     	= _newIsLocal;
				
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
				insertParticleRecursor(_newPayload, _newIsLocal);
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
				bool	residentIsLocal	= CurNodePtr->isLocal;
				
				converttoCell(targetOctant);

				insertParticleRecursor(residentPayload, residentIsLocal);
				insertParticleRecursor(_newPayload, _newIsLocal);
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
			cellData.resize(cellCounter, MSIZE);
			cellProxies.resize(cellCounter);

			for (size_t i = 0; i < cellCounter; i++) {
				( cellProxies[i] ).setup(&cellData, i);
				cellData(i, CID) = i;
			}
			std::cout << cellCounter << " cells in tree \n";
			cellCounter = 0;
			goRoot();
			connectCellDataRecursor();
			
			debugCounter = 0;
			goRoot();
			//std::cout << CurNodePtr->xCom << "   " << CurNodePtr->isEmpty << "\n";
            //std::cout << (*(CurNodePtr->payload))(CX) << "   " << CurNodePtr->isEmpty << "\n";
        
            calcMultipoleRecursor();
			
            goRoot();
            //std::cout << CurNodePtr->xCom << "   " << CurNodePtr->isEmpty << "\n";
            //std::cout << (*(CurNodePtr->payload))(CX) << "   " << CurNodePtr->isEmpty << "\n";
            
			globalSumupMultipoles();
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
            return CurNodePtr->isEmpty; // particles are per definition
                                        // empty, so we omit this check
		};
		
		/**
		* preorder recursive function to connect
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
				CurNodePtr->payload = *(cellProxies[cellCounter]);
				cellCounter++;
			}
		}

		/**
		* calculate multipole from children
		*/
		valueType monopolCM, monopolCXM, monopolCYM, monopolCZM;
		void calcMultipole() {
            monopolCM = 0.;
            monopolCXM = 0.;
            monopolCYM = 0.;
            monopolCZM = 0.;
            /**
            * check locality of current cell node:
            *  - if node is in the toptree, the node is local
            *  - if node is below toptree, the parent node is local if
            *    any child is local ( ... or ... or ... or ... )
            *
            * to check the proper working together of the costzone and
            * the toptree would be to check, that cells below toptreeDepth
            * either have only local or only non-local children, but never
            * both.
            */
            if ( ( CurNodePtr->depth ) > toptreeDepth ) {
                for (size_t i = 0; i < 8; i++) {
                    if ( CurNodePtr->child[i] != NULL ) {
                        (CurNodePtr->isLocal) |= CurNodePtr->child[i]->isLocal;
                    }
                }
            } else {
                CurNodePtr->isLocal = true;
            }
            
            /**
            * add up the contributions from the children with
            * the same locality as the current node.
            * all this locality business guarantees, that every particle
            * contributes only ONCE to the multipole moments of the global tree
            * but "ghost" cells still 
            */
            for (size_t i = 0; i < 8; i++) {
                if ( CurNodePtr->child[i] != NULL &&
                    ( CurNodePtr->child[i]->isLocal == CurNodePtr->isLocal )
                    ) {
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
										
            /* copy data to node itself ... */
            CurNodePtr->q000 = monopolCM;
            CurNodePtr->xCom = monopolCXM / monopolCM;
            CurNodePtr->yCom = monopolCYM / monopolCM;
            CurNodePtr->zCom = monopolCZM / monopolCM;
            
            /* ... and to the matrix. THE MATRIX!!! */
            (*(CurNodePtr->payload))(Q000) = CurNodePtr->q000;
            (*(CurNodePtr->payload))(CX) = CurNodePtr->xCom;
            (*(CurNodePtr->payload))(CY) = CurNodePtr->yCom;
            (*(CurNodePtr->payload))(CZ) = CurNodePtr->zCom;
		}
		
		/**
		* sum up the multipole cellData matrix globally
		* and distribute it again to each node. if MPI is not
		* defined, nothing is done here.
		*/
		void globalSumupMultipoles() {
#ifdef OOSPH_MPI
		
			matrixType recvBuffer;

			const size_t RANK = MPI::COMM_WORLD.Get_rank();
			const size_t SIZE = MPI::COMM_WORLD.Get_size();

			size_t round = 0;
			size_t remNodes = SIZE;

			const size_t noCellBytes = noToptreeCells*MSIZE*sizeof(valueType);
            
			std::queue<size_t> sumUpSend, sumUpRecv;
			std::stack<size_t> distrSend, distrRecv;

			recvBuffer.resize( noToptreeCells, MSIZE );
            bitsetType recvCellIsFilled( noToptreeCells );
			
            goRoot();
            toptreeCounter = 0;
            cellToIsFilledVectorRecursor();
            //std::cout << cellIsFilled << "\n";
            //std::cout << cellIsFilled.count() << " toptree cells are filled " << "\n";
            //std::cout << cellIsFilled.size() << "\n";
            //std::cout << recvCellIsFilled.size() << "\n";

            
			/**
			* magic algorithm which prepares sending and receiving queues
			* for the summing up step and sending and receiving stacks for
			* the distributing step.
			*/
			while ( remNodes > 1 ) {
				size_t noPairs = lrint( floor( remNodes / 2. ) );
				remNodes -= noPairs;
				size_t stepToNext = ( 1 << round );
		
				for ( size_t i = 0; i < 2*noPairs; i += 2) {
					size_t sendRank = ( SIZE - 1 ) - stepToNext*(i+1);
					size_t recvRank = ( SIZE - 1 ) - stepToNext*i;
				
					if ( RANK == sendRank ) {
						sumUpSend.push( recvRank );
						distrRecv.push( recvRank );
					} else if ( RANK == recvRank ) {
						sumUpRecv.push( sendRank );
						distrSend.push( sendRank );
					} else {}
				}
				round++;
			}
			
			// fill cellIsFilled bitset
			
            using namespace boost::posix_time;
            ptime TimeStart, TimeStop;
            std::string timeFileName = "timing_rank";
            timeFileName += boost::lexical_cast<std::string>(RANK);
            std::fstream timeFile;
            timeFile.open(timeFileName.c_str(), std::ios::out);
            
            TimeStart = microsec_clock::local_time();
                                                            
			//MPI::COMM_WORLD.Barrier();
			/**
			* receive multipoles from other nodes and add
			* them to local value
			*/
			while ( ! sumUpRecv.empty() ) {                
				size_t recvFrom = sumUpRecv.front();
				sumUpRecv.pop();

                timeFile  << "  RANK " << RANK
                          << " <- " << recvFrom << "  "
                          << " recv  bitset        "
                          << ( microsec_clock::local_time() - TimeStart )
                          << "\n";
                
                recvBitset( recvCellIsFilled, recvFrom );

                timeFile  << "  RANK " << RANK
                          << " <- " << recvFrom << "  "
                          << " recvd bitset        "
                          << ( microsec_clock::local_time() - TimeStart )
                          << "\n";
                
				MPI::COMM_WORLD.Recv( &recvBuffer(0, 0) , noCellBytes,
					MPI_BYTE, recvFrom, RANK );
                    
                timeFile  << "  RANK " << RANK
                          << " <- " << recvFrom << "  "
                          << " recvd multipoles    "
                          << ( microsec_clock::local_time() - TimeStart )
                          << "\n";
            
			// add buffer to local value
			// cellData += recvBuffer
			}
			
			/**
			* send local value to another node
			*/
			while ( ! sumUpSend.empty() ) {
				size_t sendTo = sumUpSend.front();
				sumUpSend.pop();

                timeFile  << "  RANK " << RANK
                          << " -> " << sendTo << "  "
                          << " send bitset         "
                          << ( microsec_clock::local_time() - TimeStart )
                          << "\n";
                
                sendBitset( cellIsFilled, sendTo );

                timeFile  << "  RANK " << RANK
                          << " -> " << sendTo << "  "
                          << " sent bitset         "
                          << ( microsec_clock::local_time() - TimeStart )
                          << "\n";

				/**
				 * do a synchronous send, which is non-blocking but
				 * still prevents the receiving node from getting
				 * bombarded by multiple sends
				 * */
				MPI::COMM_WORLD.Ssend( &cellData(0, 0), noCellBytes,
					MPI_BYTE, sendTo, sendTo );

                timeFile  << "  RANK " << RANK
                          << " -> " << sendTo << "  "
                          << " sent multipoles     "
                          << ( microsec_clock::local_time() - TimeStart )
                          << "\n";

			}
			
			/**
			* receive global result
			*/
			while ( ! distrRecv.empty() ) {
				size_t recvFrom = distrRecv.top();
				distrRecv.pop();

                timeFile  << "  RANK " << RANK
                          << " <- " << recvFrom << "  "
                          << " wait f. gl. multip. "
                          << ( microsec_clock::local_time() - TimeStart )
                          << "\n";


				MPI::COMM_WORLD.Recv( &cellData(0, 0), noCellBytes,
					MPI_BYTE, recvFrom, RANK );

                timeFile  << "  RANK " << RANK
                          << " <- " << recvFrom << "  "
                          << " recvd. gl. multip.  "
                          << ( microsec_clock::local_time() - TimeStart )
                          << "\n";

            }
			
			/**
			* ... and distribute it to other nodes
			*/
			while ( ! distrSend.empty() ) {
				size_t sendTo = distrSend.top();
				distrSend.pop();

				/**
				 * do a synchronous send, which is non-blocking but
				 * still prevents the receiving node from getting
				 * bombarded by multiple sends
				 * */
                timeFile  << "  RANK " << RANK
                          << " -> " << sendTo << "  "
                          << " sent global bitset  "
                          << ( microsec_clock::local_time() - TimeStart )
                          << "\n";

				MPI::COMM_WORLD.Ssend( &cellData(0, 0), noCellBytes,
					MPI_BYTE, sendTo, sendTo );
                
                timeFile  << "  RANK " << RANK
                          << " -> " << sendTo << "  "
                          << " sent global multip. "
                          << ( microsec_clock::local_time() - TimeStart )
                          << "\n";

			}
#endif		
		}
		
		/**
		* preorder recursive function to connect
		* cells to a matrix row
		*/
        
        size_t toptreeCounter;
		void cellToIsFilledVectorRecursor(void) {
            if ( belowToptreeStop() ) {} else {
                cellIsFilled[toptreeCounter] = !(CurNodePtr->isEmpty);
                /*std::cout << cellIsFilled[toptreeCounter] << "   "
                          << CurNodePtr->q000 << "   "
                          << CurNodePtr->depth << "\n";*/
                toptreeCounter++;                        
            
                for (size_t i = 0; i < 8; i++) { // try without loop
                    if ( CurNodePtr->child[i] != NULL ) {
                        goChild(i);
                        
                        cellToIsFilledVectorRecursor();
                        goUp();
                    }
                }
            }
		}

		void isFilledVectorToCellRecursor(void) {
            if ( belowToptreeStop() ) {} else {
                CurNodePtr->isEmpty = ! cellIsFilled[toptreeCounter];
                toptreeCounter++;
                                
                for (size_t i = 0; i < 8; i++) { // try without loop
                    if ( CurNodePtr->child[i] != NULL ) {
                        goChild(i);
                        isFilledVectorToCellRecursor();
                        goUp();
                    }
                }
            }
		}

		/**
		* stop recursion if:
		* - depth of current node is below toptreeDepth
		*/
		bool belowToptreeStop(void) {
			return ( CurNodePtr->depth > toptreeDepth );
		};
        
        /**
        * little comm helper
        * \todo move to comm_manager later on
        */
        
        void recvBitset( bitsetRefType _bitSet, size_t _sendRank ) {
            // prepare buffer
            const size_t RANK = MPI::COMM_WORLD.Get_rank();
            size_t noBlocks = lrint( ceil(
                (double)noToptreeCells
                / ( sizeof(bitsetBlockType) * 8 ) ) );
            size_t noBsBytes = noBlocks*sizeof(bitsetBlockType);
            std::vector<bitsetBlockType> recvBuff(noBlocks);
            
            // receive
            MPI::COMM_WORLD.Recv( &(recvBuff[0]), noBsBytes,
                MPI_BYTE, _sendRank, RANK );
                
            // copy back to bitset
            boost::from_block_range( recvBuff.begin(), recvBuff.end(),
                _bitSet );
        }
        
        void sendBitset( bitsetRefType _bitSet, size_t _recvRank ) {
            // prepare buffer
            const size_t RANK = MPI::COMM_WORLD.Get_rank();
            size_t noBlocks = lrint( ceil(
                (double)noToptreeCells
                / ( sizeof(bitsetBlockType) * 8 ) ) );
            size_t noBsBytes = noBlocks*sizeof(bitsetBlockType);
            std::vector<bitsetBlockType> sendBuff(noBlocks);
            
            // copy bitset to buffer
            boost::to_block_range( _bitSet, sendBuff.begin() );
            
            std::cout << sendBuff.size() << " " << _bitSet.size() << " " << noBsBytes << "\n";
            
            // send
            MPI::COMM_WORLD.Ssend( &(sendBuff[0]), noBsBytes,
                MPI_BYTE, _recvRank, _recvRank );
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

