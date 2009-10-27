#ifndef BHTREE_NODES_CPP
#define BHTREE_NODES_CPP

/*
 *  bhtree_nodes.cpp
 *
 *
 *  Created by Andreas Reufer on 07.10.09
 *  Copyright 2009 University of Berne. All rights reserved.
 *
 */

#include "bhtree_nodes.h"
#include "bhtree_particle.h"

namespace sphlatch {
///
/// pointer operator for nodes
///
nodePtrT genericNode::operator*()
{
   return(this);
}

///
/// clear methods for the various cell node types
///
void genericNode::clear()
{
   parent = NULL;
   next   = NULL;

   depth = -1;
   ident = 0;

   isParticle = false;
   isCZ       = false;
   isRemote   = false;

   atBottom = false;
   neighSet = false;

   isSettled   = false;
   needsUpdate = true;
}

void genericCellNode::clear()
{
   genericNode::clear();

   skip = NULL;
   for (size_t i = 0; i < 8; i++)
   {
      child[i] = NULL;
   }

   cen  = 0., 0., 0.;
   clSz = 0.;
}

void monopoleCellNode::clear()
{
   genericCellNode::clear();

   com  = 0., 0., 0.;
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

   domain  = 0;
   absCost = 0.;
   noParts = 0;

   orphFrst = NULL;
   orphLast = NULL;
   chldFrst = NULL;
   chldLast = NULL;

   isCZ = true;
}

void particleNode::clear()
{
   genericNode::clear();

   isParticle = true;

   partPtr = NULL;
   pos     = 0., 0., 0.;
   m       = 0.;
}

///
/// generic cell methods
///
void genericCellNode::inheritCellPos(size_t _n)
{
   clSz = 0.5 * (static_cast<gcllPtrT>(parent)->clSz);

   const fType hcSize = 0.5 * clSz;

   cen = static_cast<gcllPtrT>(parent)->cen;

   cen[0] += ((_n >> 0) % 2) ? hcSize : -hcSize;
   cen[1] += ((_n >> 1) % 2) ? hcSize : -hcSize;
   cen[2] += ((_n >> 2) % 2) ? hcSize : -hcSize;

   depth = parent->depth + 1;
}

bool genericCellNode::pointInsideCell(const vect3dT& _pos)
{
   const fType hclSz = 0.5 * clSz;

   //return ( all( cen - _pos < hclSz ) && all( cen - _pos > -hclSz ) );
   return(cen[0] - hclSz < _pos[0] && cen[0] + hclSz > _pos[0] &&
          cen[1] - hclSz < _pos[1] && cen[1] + hclSz > _pos[1] &&
          cen[2] - hclSz < _pos[2] && cen[2] + hclSz > _pos[2]);
}

size_t genericCellNode::getOctant(const vect3dT& _pos)
{
   size_t targetOctant = 0;

   targetOctant += _pos[0] < cen[0] ? 0 : 1;
   targetOctant += _pos[1] < cen[1] ? 0 : 2;
   targetOctant += _pos[2] < cen[2] ? 0 : 4;
   return(targetOctant);
}

///
/// methods to init a costzone from a quadrupole cell and vice versa
///
void quadrupoleCellNode::initFromCZll(czllT& _czll)
{
   *this    = static_cast<qcllT>(_czll);
   isCZ     = false;
   atBottom = false;
   neighSet = false;

   for (size_t i = 0; i < 8; i++)
   {
      if (child[i] != NULL)
         child[i]->parent = this;
   }

   isSettled   = false;
   needsUpdate = true;
}

///
/// calculate the quadrupole moments of a cell
///
void quadrupoleCellNode::calcMultipole()
{ }

void costzoneCellNode::initFromCell(qcllT& _cell)
{
   *static_cast<qcllPtrT>(this) = _cell;
   isCZ = true;

   for (size_t i = 0; i < 8; i++)
   {
      if (child[i] != NULL)
         child[i]->parent = this;
   }

   for (size_t i = 0; i < 27; i++)
   {
      neighbour[i] = NULL;
   }

   domain   = 0;
   atBottom = false;
   neighSet = false;

   absCost = 0.;
}

///
/// push down the neighbours to the childs
//FIXME: this is heavily untested!
///
void costzoneCellNode::pushdownNeighbours()
{
   ///
   /// if the current cells neighbours are not set,
   /// then push down neighbours of parent
   ///
   if (not neighSet)
      static_cast<czllPtrT>(parent)->pushdownNeighbours();

   ///
   /// temporary neighbour pointer
   ///
   nodePtrT neighPtr = NULL;

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

                        static_cast<czllPtrT>(child[jc])->neighbour[js] =
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
                  if ((static_cast<gcllPtrT>(neighbour[ps])->clSz == clSz) &&
                      not static_cast<czllPtrT>(neighbour[ps])->atBottom)
                  {
                     ///
                     /// determine the neighbours child index
                     ///
                     const int jc = ((ix + 2) % 2) +
                                    ((iy + 2) % 2) * 2 +
                                    ((iz + 2) % 2) * 4;

                     if (static_cast<gcllPtrT>(neighbour[ps])->child[jc]
                         != NULL)
                        neighPtr =
                           static_cast<gcllPtrT>(neighbour[ps])->child[jc];
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

                        static_cast<czllPtrT>(child[jc])->neighbour[js] =
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

///
/// adopt a particle
///
void costzoneCellNode::adopt(pnodPtrT _pnod)
{
   if (orphFrst == NULL)
      orphFrst = _pnod;
   else
      orphLast->next = _pnod;

   _pnod->next = NULL;
   orphLast    = _pnod;
}

///
/// update a particles position
///
void particleNode::update()
{
   pos = partPtr->pos;
   m   = partPtr->m;
}

///
/// deletes a CZ cells subtree
/// (the childs walk has to be set!)
///
void costzoneCellNode::delSubtree()
{
   nodePtrT curChld = chldFrst, nxtChld = NULL;

   chldLast->next = NULL;

   while (curChld != NULL)
   {
      nxtChld = curChld->next;
      if (curChld->isParticle)
         delete static_cast<pnodPtrT>(curChld);
      else if (curChld->isCZ)
         delete static_cast<czllPtrT>(curChld);
      else
         delete static_cast<qcllPtrT>(curChld);
      curChld = nxtChld;
   }
}
};
#endif
