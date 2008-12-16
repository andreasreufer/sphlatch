#ifndef BHTREE_NODES_H
#define BHTREE_NODES_H

//#define SPHLATCH_AMD64PADDING

/*
 *  bhtree_nodes.h
 *
 *
 *  Created by Andreas Reufer on 26.11.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */


#include "typedefs.h"

namespace sphlatch {
///
/// base class for all nodes
///
class genericNode {
public:
   genericNode* parent;
   genericNode* next;
   genericNode* skip;
   size_t     depth;

   identType ident;

   bool isParticle : 1;
   bool isRemote   : 1;
   bool isCZ       : 1;

   genericNode* operator*()
   {
      return(this);
   }

   genericNode()
   {
      clear();
   }

   ~genericNode() { }

   void clear()
   {
      parent = NULL;
      next   = NULL;
      skip   = NULL;

      depth = -1;
      ident = 0;

      isParticle = false;
      isRemote   = false;
      isCZ       = false;
   }
};

///
/// particle node
///
class particleNode : public genericNode {
public:

   particleNode()
   {
      clear();
   }

   void clear()
   {
      parent = NULL;
      next   = NULL;
      skip   = NULL;

      depth = -1;
      ident = 0;

      isParticle = true;
      isRemote   = false;
      isCZ       = false;

      xPos = 0.;
      yPos = 0.;
      zPos = 0.;
      mass = 0.;
   }

   fType xPos, yPos, zPos;
   fType mass;
#ifdef SPHLATCH_AMD64PADDING
private:
   char pad[0];
#endif
};

///
/// generic cell node
///
class genericCellNode : public genericNode {
public:
   typedef genericNode*       genericNodePtrT;
   typedef genericCellNode*   selfPtr;

   genericNodePtrT child[8];
   fType           xCen, yCen, zCen;
   fType           clSz;
   fType           cost;

   ///
   /// generic cell node constructor
   ///
   genericCellNode()
   {
      clear();
   }

   ~genericCellNode()
   {
   }

   void clear()
   {
      parent = NULL;
      next   = NULL;
      skip   = NULL;

      depth = -1;
      ident = 0;

      isParticle = false;
      isRemote   = false;
      isCZ       = false;

      for (size_t i = 0; i < 8; i++)
      {
         child[i] = NULL;
      }

      xCen = 0.;
      yCen = 0.;
      zCen = 0.;
      clSz = 0.;
      cost = 0.;
   }

   void setCellPos(size_t _n)
   {
      clSz = 0.5 * (static_cast<selfPtr>(parent)->clSz);
      const fType hcSize = 0.5 * clSz;

      xCen = ((_n) % 2) ?
             static_cast<selfPtr>(parent)->xCen + hcSize :
             static_cast<selfPtr>(parent)->xCen - hcSize;
      yCen = ((_n >> 1) % 2) ?
             static_cast<selfPtr>(parent)->yCen + hcSize :
             static_cast<selfPtr>(parent)->yCen - hcSize;
      zCen = ((_n >> 2) % 2) ?
             static_cast<selfPtr>(parent)->zCen + hcSize :
             static_cast<selfPtr>(parent)->zCen - hcSize;
   }
};

///
/// cell node with monopole moments
///
class monopoleCellNode : public genericCellNode {
public:
   monopoleCellNode()
   {
      clear();
   }

   ~monopoleCellNode() { }

   void clear()
   {
      parent = NULL;
      next   = NULL;
      skip   = NULL;

      depth = -1;
      ident = 0;

      isParticle = false;
      isRemote   = false;
      isCZ       = false;

      for (size_t i = 0; i < 8; i++)
      {
         child[i] = NULL;
      }

      xCen = 0.;
      yCen = 0.;
      zCen = 0.;
      clSz = 0.;
      cost = 0.;

      xCom = 0.;
      yCom = 0.;
      zCom = 0.;
      mass = 0.;
   }

   fType xCom, yCom, zCom;
   fType mass;
};

///
/// cell node with quadrupole moments
///
class quadrupoleCellNode : public monopoleCellNode {
public:
   fType q11, q22, q33, q12, q13, q23;

   quadrupoleCellNode()
   {
      clear();
   }

   ~quadrupoleCellNode() { }

   void clear()
   {
      parent = NULL;
      next   = NULL;
      skip   = NULL;

      depth = -1;
      ident = 0;

      isParticle = false;
      isRemote   = false;
      isCZ       = false;

      for (size_t i = 0; i < 8; i++)
      {
         child[i] = NULL;
      }

      xCen = 0.;
      yCen = 0.;
      zCen = 0.;
      clSz = 0.;
      cost = 0.;

      xCom = 0.;
      yCom = 0.;
      zCom = 0.;
      mass = 0.;

      q11 = 0.;
      q22 = 0.;
      q33 = 0.;
      q12 = 0.;
      q13 = 0.;
      q23 = 0.;
   }

   void copy(const quadrupoleCellNode* _cellPtr)
   {
      parent = _cellPtr->parent;
      next   = _cellPtr->next;
      skip   = _cellPtr->skip;

      depth = _cellPtr->depth;
      ident = _cellPtr->ident;

      isParticle = false;
      isRemote   = _cellPtr->isRemote;
      isCZ       = false;

      for (size_t i = 0; i < 8; i++)
      {
         child[i] = _cellPtr->child[i];
      }

      xCen = _cellPtr->xCen;
      yCen = _cellPtr->yCen;
      zCen = _cellPtr->zCen;
      clSz = _cellPtr->clSz;
      cost = _cellPtr->cost;

      xCom = _cellPtr->xCom;
      yCom = _cellPtr->yCom;
      zCom = _cellPtr->zCom;
      mass = _cellPtr->mass;

      q11 = _cellPtr->q11;
      q22 = _cellPtr->q22;
      q33 = _cellPtr->q33;
      q12 = _cellPtr->q12;
      q13 = _cellPtr->q13;
      q23 = _cellPtr->q23;
   }

private:
#ifdef SPHLATCH_AMD64PADDING
   char pad[48];
#endif
};

///
/// a quadrupole costzone cell
///
class costzoneCellNode : public quadrupoleCellNode {
public:
   typedef costzoneCellNode*   selfPtr;

   genericNodePtrT sibling[27];
   identType       domain;
   bool            atBottom : 1;

   costzoneCellNode()
   {
      clear();
   }

   ~costzoneCellNode() { }

   void clear()
   {
      parent = NULL;
      next   = NULL;
      skip   = NULL;

      depth = -1;
      ident = 0;

      isParticle = false;
      isRemote   = false;
      isCZ       = true;

      for (size_t i = 0; i < 8; i++)
      {
         child[i] = NULL;
      }

      xCen = 0.;
      yCen = 0.;
      zCen = 0.;
      clSz = 0.;
      cost = 0.;

      xCom = 0.;
      yCom = 0.;
      zCom = 0.;
      mass = 0.;

      q11 = 0.;
      q22 = 0.;
      q33 = 0.;
      q12 = 0.;
      q13 = 0.;
      q23 = 0.;

      for (size_t i = 0; i < 27; i++)
      {
         sibling[i] = NULL;
      }

      domain   = 0;
      atBottom = false;
   }

   ///
   /// inherit siblings from parents and try to push
   /// them down the tree
   ///
   void inheritSiblings()
   {
      domain = static_cast<selfPtr>(parent)->domain;
      for (size_t i = 0; i < 27; i++)
      {
         sibling[i] = static_cast<selfPtr>(parent)->sibling[i];

         ///
         /// check whether we can pull it down
         ///
      }
      sibling[13] = this;
   }

private:
#ifdef SPHLATCH_AMD64PADDING
   char pad[28];
#endif
};
};
#endif
