#ifndef BHTREE_PARTICLE_H
#define BHTREE_PARTICLE_H

/*
 *  bhtree_particle.h
 *
 *
 *  Created by Andreas Reufer on 23.02.09
 *  Copyright 2009 University of Berne. All rights reserved.
 *
 */

#include "bhtree_nodes.h"
#include "typedefs.h"

namespace sphlatch {
class particleNode;
///
/// basic tree particle class
///
class treeGhost {
public:
   typedef particleNode*   partNodePtrT;
   vect3dT      pos;
   fType        m, eps;
   partNodePtrT treeNode;

   idType id;
   fType  cost;

#ifdef SPHLATCH_PADD64
private:
   char pad[8];
#endif
};

///
/// dummy resident class
///
class treePart : public treeGhost { };

///
/// moving ghost particle class
///
class movingGhost {
  public:
   vect3dT vel;
};

///
/// moving resident particle class
///
class movingPart : public movingGhost {
  public:
   vect3dT acc;
};
};
#endif
