#ifndef PARTICLE_H
#define PARTICLE_H

//#define SPHLATCH_PADTO64BYTES

/*
 *  particle.h
 *
 *
 *  Created by Andreas Reufer on 23.02.09
 *  Copyright 2009 University of Berne. All rights reserved.
 *
 */


#include "typedefs.h"

namespace sphlatch {
class particleNode;

class treeParticle {
public:
   typedef particleNode*   partNodePtrT;
   vect3dT      pos;
   fType        m, eps;
   partNodePtrT treeNode;

   idType id;
   fType  cost;

#ifdef SPHLATCH_PADTO64BYTES
private:
   char pad[8];
#endif
};

class SPHghost {
public:
   fType h, rho;
#ifdef SPHLATCH_PADTO64BYTES
private:
   char pad[8];
#endif
};

class SPHpart : public SPHghost {
public:
  fType dhdt;
#ifdef SPHLATCH_PADTO64BYTES
private:
   char pad[0];
#endif
};

};
#endif
