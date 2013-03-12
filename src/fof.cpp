#ifndef SPHLATCH_FOF_CPP
#define SPHLATCH_FOF_CPP

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
#include "friend_particle.h"
#include "clump.cpp"

#include "bhtree.cpp"
#include "bhtree_worker_neighfunc.cpp"

namespace sphlatch {
template<typename _partT>
class FOF : public ParticleSet<Clump<_partT> >
{
public:
   typedef _partT                 partT;
   typedef ParticleSet<_partT>    partSetT;
   typedef std::list<_partT*>     partPtrLT;

   typedef BHTree                 treeT;

   typedef Clump<_partT>          friendT;
   typedef ParticleSet<friendT>   parentT;

   typedef std::list<friendT>     clumpLT;


public:
   FOF(partSetT& _parts,
       treeT&    _tree) : parts(_parts), tree(_tree),
                          friends(parentT::parts)
   {
      nextFreeId          = 1;
      parentT::step       = _parts.step;
      parentT::attributes = _parts.attributes;
   }

   ~FOF() { }

   void search(const fType _rhoMin, const fType _hMult, const fType _mMin);

private:
   cType mergeFOF(const cType cida, const cType cidb);

   partSetT& parts;
   treeT&    tree;

   std::vector<friendT>& friends;

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
      bool operator()(const friendT& _c1, const friendT& _c2)
      {
         return(_c1.m > _c2.m);
      }
   };

   cType nextFreeId;
};

template<typename _partT>
cType FOF<_partT>::mergeFOF(const cType cida, const cType cidb)
{
   if (cida == cidb)
      return(cida);

   cType nid, oid;

   if (cida < cidb)
   {
      nid = cida;
      oid = cidb;
   }
   else
   {
      nid = cidb;
      oid = cida;
   }

   const size_t nop = parts.getNop();
   for (size_t i = 0; i < nop; i++)
   {
      if (parts[i].friendid == oid)
         parts[i].friendid = nid;
   }

   return(nid);
}

template<typename _partT>
void FOF<_partT>::search(const fType _rhoMin, const fType _hMult,
                         const fType _mMin)
{
   const size_t nop = parts.getNop();

   for (size_t i = 0; i < nop; i++)
      parts[i].friendid = FRIENDNOTSET;

   NeighFindWorker<_partT> NFW(&tree);

   for (size_t i = 0; i < nop; i++)
   {
      partT* cPart = &(parts[i]);

      //
      // if density is too low, particle does not belong to any clump
      // regardless of previous state
      //
      if (parts[i].rho < _rhoMin)
      {
         parts[i].friendid = FRIENDNONE;
         continue;
      }
      else
      {
         if ((parts[i].friendid == FRIENDNONE) ||
             (parts[i].friendid == FRIENDNOTSET))
         {
            const fType srad   = _hMult * cPart->h;
            partPtrLT   neighs = NFW(cPart, srad);

            std::set<cType> neighFOF;
            typename partPtrLT::iterator nitr;
            for (nitr = neighs.begin(); nitr != neighs.end(); nitr++)
               if (((*nitr)->rho > _rhoMin) && ((*nitr)->friendid > 0))
                  neighFOF.insert((*nitr)->friendid);

            const size_t noFOFFound = neighFOF.size();
            switch (noFOFFound)
            {
            case 0:
            {
               // so all particles
               cType nfriendid = nextFreeId;
               nextFreeId++;

               partPtrLT clumpees;
               for (nitr = neighs.begin(); nitr != neighs.end(); nitr++)
                  if ((*nitr)->rho > _rhoMin)
                  {
                     (*nitr)->friendid = FRIENDONLIST;
                     clumpees.push_back(*nitr);
                  }
                  else
                     (*nitr)->friendid = FRIENDNONE;

               // start a new clump.
               size_t nocp = 0;
               while (not clumpees.empty())
               {
                  partT* curf = clumpees.back();
                  clumpees.pop_back();

                  curf->friendid = nfriendid;
                  nocp++;

                  neighs.clear();
                  neighs = NFW(curf, curf->h);

                  for (nitr = neighs.begin(); nitr != neighs.end(); nitr++)
                  {
                     if ((*nitr)->rho > _rhoMin)
                     {
                        if ((*nitr)->friendid == nfriendid)
                           continue;

                        switch ((*nitr)->friendid)
                        {
                        case FRIENDNONE:
                           (*nitr)->friendid = FRIENDONLIST;
                           clumpees.push_back(*nitr);
                           break;

                        case FRIENDNOTSET:
                           (*nitr)->friendid = FRIENDONLIST;
                           clumpees.push_back(*nitr);
                           break;

                        case FRIENDONLIST:
                           break;

                        default:
                           nfriendid = mergeFOF((*nitr)->friendid,
                                                nfriendid);
                           (*nitr)->friendid = nfriendid;
                           break;
                        }
                     }
                     else
                     {
                        (*nitr)->friendid = FRIENDNONE;
                     }
                  }
               }
            }
               break;

            case 1:
            {
               cPart->friendid = *(neighFOF.begin());
            }
               break;

            default:
            {
               cType nfriendid = *(neighFOF.begin());

               std::set<cType>::const_iterator ncItr;
               for (ncItr = neighFOF.begin(); ncItr != neighFOF.end();
                    ncItr++)
                  nfriendid = mergeFOF(*ncItr, nfriendid);

               cPart->friendid = nfriendid;
            }
               break;
            }
         }
         else
         { }
      }
   }

   // now collect the friends
   const size_t noc = nextFreeId;

   std::vector<friendT> massRanked(noc);
   for (size_t i = 0; i < noc; i++)
      massRanked[i].id = i;

   for (size_t i = 0; i < nop; i++)
   {
      if (parts[i].friendid < FRIENDNONE)
         parts[i].friendid = FRIENDNONE;

      if (parts[i].friendid > FRIENDNONE)
      {
         const size_t cID = parts[i].friendid;
         assert(cID < noc);
         massRanked[cID].addParticle(&parts[i]);
      }
   }

   std::sort(massRanked.begin(), massRanked.end(), higherClumpMass());

   std::vector<size_t> massRankedInv(noc);
   for (size_t i = 0; i < noc; i++)
      massRankedInv[massRanked[i].id] = i;

   size_t maxID = 0.;
   for (size_t i = 0; i < noc; i++)
   {
      if (massRanked[i].m < _mMin)
         break;
      maxID = i + 1;
   }

   for (size_t i = 0; i < nop; i++)
   {
      const size_t oldID = parts[i].friendid;
      const size_t newID = massRankedInv[oldID] + 1;

      if (newID <= maxID)
         parts[i].friendid = newID;
      else
         parts[i].friendid = FRIENDNONE;
   }

   const size_t fnoc = maxID + 1;
   parentT::resize(fnoc);

   for (size_t i = 0; i < nop; i++)
      friends[parts[i].friendid].addParticle(&(parts[i]));


   for (size_t i = 0; i < fnoc; i++)
   {
      friendT& curClump(parentT::operator[](i));
      curClump.calcProperties();
      curClump.calcAngularMom();
      curClump.assignID(i);
   }
   
   return;
}
}

#endif
