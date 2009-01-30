#ifndef BHTREE_NODE_CELLS_H
#define BHTREE_NODE_CELLS_H

/*
 *  bhtree_node_cells.h
 *
 *
 *  Created by Andreas Reufer on 20.01.09
 *  Copyright 2009 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"
#include "bhtree_node_generic.h"

namespace sphlatch {
///
/// generic cell node
///
class genericCellNode : public genericNode {
public:
   typedef genericNode*       nodePtr;
   typedef genericCellNode*   gcllPtr;

   nodePtr child[8];

   fType xCen, yCen, zCen;
   fType clSz, cost;

   genericCellNode() { }
   ~genericCellNode() { }

   void clear();
   void inheritCellPos(size_t _n);
   bool pointInsideCell(const fType _x, const fType _y, const fType _z);
   size_t getOctant(const fType _x, const fType _y, const fType _z);

private:
};

///
/// cell node with monopole moments
///
class monopoleCellNode : public genericCellNode {
public:
   fType xCom, yCom, zCom;
   fType mass;

   monopoleCellNode() { }
   ~monopoleCellNode() { }

   void clear();

#ifdef SPHLATCH_AMD64PADDING
   char pad[0];
#endif
};

///
/// cell node with quadrupole moments
///
class costzoneCellNode;
class quadrupoleCellNode : public monopoleCellNode {
public:
   typedef costzoneCellNode     czllT;
   typedef quadrupoleCellNode   cellT;

   fType q11, q22, q33, q12, q13, q23;

   quadrupoleCellNode() { }
   ~quadrupoleCellNode() { }

   void clear();
   void initFromCZll(czllT& _czll);

private:
#ifdef SPHLATCH_AMD64PADDING
   char pad[16];
#endif
};

///
/// a quadrupole costzone cell
///
class costzoneCellNode : public quadrupoleCellNode {
public:
   typedef costzoneCellNode*     czllPtr;
   typedef quadrupoleCellNode*   cellPtr;

   genericNodePtrT neighbour[27];
   identType       domain;
   countsType       noParts;

   costzoneCellNode() { }
   ~costzoneCellNode() { }

   void clear();
   void initFromCell(cellT& _cell);

   void pushdownNeighbours();

private:
#ifdef SPHLATCH_AMD64PADDING
   char pad[12];
#endif
};


///
/// clear functions for the various cell node types
///
void genericCellNode::clear()
{
   genericNode::clear();

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

void monopoleCellNode::clear()
{
   genericCellNode::clear();

   xCom = 0.;
   yCom = 0.;
   zCom = 0.;
   mass = 0.;
}

void quadrupoleCellNode::clear()
{
   monopoleCellNode::clear();

   q11 = 0.;
   q22 = 0.;
   q33 = 0.;
   q12 = 0.;
   q13 = 0.;
   q23 = 0.;
}

void costzoneCellNode::clear()
{
   quadrupoleCellNode::clear();

   for (size_t i = 0; i < 27; i++)
   {
      neighbour[i] = NULL;
   }

   domain   = 0;
   noParts = 0;
   atBottom = false;
}

void genericCellNode::inheritCellPos(size_t _n)
{
   clSz = 0.5 * (static_cast<gcllPtr>(parent)->clSz);
   const fType hcSize = 0.5 * clSz;

   xCen = ((_n) % 2) ?
          static_cast<gcllPtr>(parent)->xCen + hcSize :
          static_cast<gcllPtr>(parent)->xCen - hcSize;
   yCen = ((_n >> 1) % 2) ?
          static_cast<gcllPtr>(parent)->yCen + hcSize :
          static_cast<gcllPtr>(parent)->yCen - hcSize;
   zCen = ((_n >> 2) % 2) ?
          static_cast<gcllPtr>(parent)->zCen + hcSize :
          static_cast<gcllPtr>(parent)->zCen - hcSize;

   depth = parent->depth + 1;
}

bool genericCellNode::pointInsideCell(const fType _x,
                                      const fType _y,
                                      const fType _z)
{
   const fType hclSz = 0.5 * clSz;

   return(xCen - hclSz < _x && xCen + hclSz > _x &&
          yCen - hclSz < _y && yCen + hclSz > _y &&
          zCen - hclSz < _z && zCen + hclSz > _z);
}

size_t genericCellNode::getOctant(const fType _x,
                                  const fType _y,
                                  const fType _z)
{
   size_t targetOctant = 0;

   targetOctant += _x < xCen ? 0 : 1;
   targetOctant += _y < yCen ? 0 : 2;
   targetOctant += _z < zCen ? 0 : 4;
   return(targetOctant);
}

///
/// methods to init a costzone from a quadrupole cell and vice versa
///
void quadrupoleCellNode::initFromCZll(czllT& _czll)
{
   *this = static_cast<cellT>(_czll);
   isCZ  = false;
}

void costzoneCellNode::initFromCell(quadrupoleCellNode& _cell)
{
   *static_cast<cellPtr>(this) = _cell;
   isCZ = true;

   for (size_t i = 0; i < 27; i++)
   {
      neighbour[i] = NULL;
   }

   domain   = 0;
   atBottom = false;
}

///
/// push down the neighbours to the childs
///
void costzoneCellNode::pushdownNeighbours()
{
   ///
   /// if the current cells neighbours are not set,
   /// then push down neighbours of parent
   ///
   if (not neighSet)
      static_cast<czllPtr>(parent)->pushdownNeighbours();

   ///
   /// temporary neighbour pointer
   ///
   nodePtr neighPtr = NULL;

   ///
   /// iterate through the 4x4x4 possible childs i surrounding
   /// the 2x2x2 childs of the cell
   ///
   for (int ix = -1; ix < 3; ix++)
   {
      for (int iy = -1; iy < 3; iy++)
      {
         for (int iz = -1; iz < 3; iz++)
         {
            ///
            /// is the current cell i a sibling?
            ///
            if ((ix > -1) && (ix < 2) &&
                (iy > -1) && (iy < 2) &&
                (iz > -1) && (iz < 2))
            {
               ///
               /// determine the child i index
               ///
               const int ic = ix + 2 * iy + 4 * iz;
               //neighPtr = child[ic];

               ///
               /// iterate through all child cells j
               /// and set them as the corresponding
               /// neighbours jc. note: also the current cell
               /// is used as its own neighbour
               ///
               for (int jx = 0; jx < 2; jx++)
               {
                  for (int jy = 0; jy < 2; jy++)
                  {
                     for (int jz = 0; jz < 2; jz++)
                     {
                        ///
                        /// determine the neighbour index
                        ///
                        const int jc = jx + 2 * jy + 4 * jz;
                        const int js = (ix - jx + 1) +
                                       (iy - jy + 1) * 3 +
                                       (iz - jz + 1) * 9;

                        static_cast<czllPtr>(child[jc])->neighbour[js] =
                           child[ic];
                     }
                  }
               }
            }
            else
            {
               ///
               /// if the corresponding parent neighbour is
               /// not null, it or its child is the neighbour of
               /// the child
               ///
               const int pix = (ix + 2) / 2;
               const int piy = (iy + 2) / 2;
               const int piz = (iz + 2) / 2;

               const int ps = pix + piy * 3 + piz * 9;

               neighPtr = NULL;
               if (neighbour[ps] != NULL)
               {
                  ///
                  /// if the neighbour cell has the same size as the
                  /// cell, it may contain a child suitable as a
                  /// neighbour. only try to use the child, when its parent
                  /// is not already at the costzone bottom.
                  ///
                  if ((static_cast<gcllPtr>(neighbour[ps])->clSz == clSz) &&
                      not static_cast<czllPtr>(neighbour[ps])->atBottom)
                  {
                     ///
                     /// determine the neighbours child index
                     ///
                     const int jc = ((ix + 2) % 2) +
                                    ((iy + 2) % 2) * 2 +
                                    ((iz + 2) % 2) * 4;

                     if (static_cast<gcllPtr>(neighbour[ps])->child[jc]
                         != NULL)
                        neighPtr =
                           static_cast<gcllPtr>(neighbour[ps])->child[jc];
                     else
                        neighPtr = neighbour[ps];
                  }
                  else
                     neighPtr = neighbour[ps];
               }

               ///
               /// iterate through all child cells j
               /// and set their corresponding neighbour
               ///
               int jxmin = std::max(0, ix - 1);
               int jymin = std::max(0, iy - 1);
               int jzmin = std::max(0, iz - 1);

               int jxmax = std::min(1, ix + 1);
               int jymax = std::min(1, iy + 1);
               int jzmax = std::min(1, iz + 1);

               for (int jx = jxmin; jx <= jxmax; jx++)
               {
                  for (int jy = jymin; jy <= jymax; jy++)
                  {
                     for (int jz = jzmin; jz <= jzmax; jz++)
                     {
                        ///
                        /// determine the neighbour index
                        ///
                        const int jc = jx + 2 * jy + 4 * jz;
                        const int js = (ix - jx + 1) +
                                       (iy - jy + 1) * 3 +
                                       (iz - jz + 1) * 9;

                        static_cast<czllPtr>(child[jc])->neighbour[js] =
                           neighPtr;
                     }
                  }
               }
            }
         }
      }
   }

   ///
   /// childs neighbours are set in this cell
   ///
   neighSet = true;
}
}

#endif
