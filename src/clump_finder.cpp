#ifndef SPHLATCH_CLUMP_FINDER_CPP
#define SPHLATCH_CLUMP_FINDER_CPP

/*
 *  clump_finder.cpp
 *
 *
 *  Created by Andreas Reufer on 10.12.09
 *  Copyright 2009 University of Berne. All rights reserved.
 *
 */

#include <list>
#include <algorithm>

#include "typedefs.h"
#include "particle_set.cpp"
#include "clump_particle.h"
#include "clump.cpp"

#include "bhtree.cpp"
#include "bhtree_worker_neighfunc.cpp"
#include "bhtree_worker_grav.cpp"

namespace sphlatch {
template<typename _partT>
class Clumps : public ParticleSet<Clump<_partT> >
{
public:
   typedef _partT                partT;
   typedef ParticleSet<_partT>   partSetT;
   typedef std::list<_partT*>    partPtrLT;

   typedef Clump<_partT>         clumpT;
   typedef ParticleSet<clumpT>   parentT;

   typedef std::list<clumpT>     clumpLT;


public:
   Clumps() : clumps(parentT::parts) { nextFreeId = 1; }
   ~Clumps() { }

   void getClumps(partSetT& _parts);

private:

   std::vector<clumpT>& clumps;

   class lowerPot
   {
public:
      bool operator()(_partT* _p1, _partT* _p2)
      {
         return(_p1->pot < _p2->pot);
      }
   };

   class higherClumpMass {
public:
      bool operator()(const clumpT& _c1, const clumpT& _c2)
      {
         return(_c1.m > _c2.m);
      }
   };

   treeT Tree;

   cType nextFreeId;
};


template<typename _partT>
void Clumps<_partT>::getClumps(ParticleSet<_partT>& _parts,
                               const fType          _minMass)
{
   const size_t        nop           = _parts.getNop();
   const fType         G             = _parts.attributes["gravconst"];
   treeT::czllPtrVectT CZbottomLoc   = Tree.getCZbottomLoc();
   const int           noCZbottomLoc = CZbottomLoc.size();

   gravT gravWorker(&Tree, G);

#pragma omp parallel for firstprivate(gravWorker)
   for (int i = 0; i < noCZbottomLoc; i++)
      gravWorker.calcAcc(CZbottomLoc[i]);

   partPtrLT unbdParts;
   for (size_t i = 0; i < nop; i++)
      unbdParts.push_back(&(_parts[i]));

   lowerPot potSorter;
   unbdParts.sort(potSorter);

   typename partPtrLT::iterator pItr, nItr, onItr;

   nextFreeId = 1;

   pItr = unbdParts.begin();

   while (pItr != unbdParts.end())
   {
      partT* cp = *pItr;

      Clump<_partT> cclmp;
      cclmp.addParticle(cp);

      size_t noFound;
      do
      {
         nItr = pItr;
         nItr++;
         if (nItr == unbdParts.end())
            break;

         noFound = 0;
         while (nItr != unbdParts.end())
         {
            if (cclmp.totEnergy(*nItr, G) < 0.)
            {
               cclmp.addParticle(*nItr);
               noFound++;
               onItr = nItr;
               nItr++;
               unbdParts.erase(onItr);
            }
            else
               nItr++;
         }
      } while (noFound > 0);

      // clump is finished now
      if (cclmp.nop > 1)
      {
         cclmp.assignID(nextFreeId);
         nextFreeId++;
         std::cout << "assigned ID " << nextFreeId << "\n";
      }


      cclmp.clear();
      pItr++;
   }
}
}

#endif
