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
   typedef _partT                        partT;
   typedef ParticleSet<_partT>           partSetT;
   typedef std::list<_partT*>            partPtrLT;

   typedef BHTree                        treeT;

   typedef Clump<_partT>                 clumpT;
   typedef ParticleSet<clumpT>           parentT;

   typedef fixThetaMAC                   macT;
   typedef GravityWorker<macT, _partT>   gravT;


public:
   Clumps() { nextFreeId = 1; }
   ~Clumps() { }

   void getClumpsFOF(partSetT& _parts, const fType _rhoMin,
                     const fType _hMult);
   void getClumpsPot(partSetT& _parts);

private:

   cType mergeClumps(const cType cida, const cType cidb, partSetT& _parts);

   void prepareSearch(partSetT& _parts);

   void collectClumps(ParticleSet<_partT>& _parts, const fType _minMass);

   void FOF(partSetT& _parts, const fType _rhoMin, const fType _hMult);
   void potSearch(partSetT& _parts);
   void potSearchOld(partSetT& _parts);

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
cType Clumps<_partT>::mergeClumps(const cType cida, const cType cidb,
                                  partSetT& _parts)
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

   const size_t nop = _parts.getNop();
   for (size_t i = 0; i < nop; i++)
   {
      if (_parts[i].clumpid == oid)
         _parts[i].clumpid = nid;
   }

   return(nid);
}

template<typename _partT>
void Clumps<_partT>::getClumpsPot(ParticleSet<_partT>& _parts)
{
   prepareSearch(_parts);
   
   potSearch(_parts);

   const fType minMass = 1.1926e+23;
   collectClumps(_parts, minMass);
}

template<typename _partT>
void Clumps<_partT>::getClumpsFOF(ParticleSet<_partT>& _parts,
                                  const fType          _rhoMin,
                                  const fType          _hMult)
{
   prepareSearch(_parts);
   
   FOF(_parts, _rhoMin, _hMult);
   
   const fType minMass = 1.1926e+23;
   collectClumps(_parts, minMass);
}

template<typename _partT>
void Clumps<_partT>::prepareSearch(ParticleSet<_partT>& _parts)
{
   parentT::step = _parts.step;

   const size_t nop       = _parts.getNop();
   const fType  costppart = 1. / nop;

   Tree.setExtent(_parts.getBox() * 1.1);

   for (size_t i = 0; i < nop; i++)
   {
      _parts[i].cost = costppart;
      Tree.insertPart(_parts[i]);
   }
   Tree.update(0.8, 1.2);

   for (size_t i = 0; i < nop; i++)
      _parts[i].clumpid = CLUMPNOTSET;

   nextFreeId = 1;
}

template<typename _partT>
void Clumps<_partT>::FOF(ParticleSet<_partT>& _parts,
                         const fType          _rhoMin,
                         const fType          _hMult)
{
   const size_t nop = _parts.getNop();

   NeighFindWorker<_partT> NFW(&Tree);
   for (size_t i = 0; i < nop; i++)
   {
      partT* cPart = &(_parts[i]);

      //
      // if density is too low, particle does not belong to any clump
      // regardless of previous state
      //
      if (_parts[i].rho < _rhoMin)
      {
         _parts[i].clumpid = CLUMPNONE;
         continue;
      }
      else
      {
         if ((_parts[i].clumpid == CLUMPNONE) ||
             (_parts[i].clumpid == CLUMPNOTSET))
         {
            const fType srad   = _hMult * cPart->h;
            partPtrLT   neighs = NFW(cPart, srad);

            std::set<cType> neighClumps;
            typename partPtrLT::iterator nitr;
            for (nitr = neighs.begin(); nitr != neighs.end(); nitr++)
               if (((*nitr)->rho > _rhoMin) && ((*nitr)->clumpid > 0))
                  neighClumps.insert((*nitr)->clumpid);

            const size_t noClumpsFound = neighClumps.size();
            switch (noClumpsFound)
            {
            case 0:
            {
               // so all particles
               cType nclumpid = nextFreeId;
               nextFreeId++;

               partPtrLT clumpees;
               for (nitr = neighs.begin(); nitr != neighs.end(); nitr++)
                  if ((*nitr)->rho > _rhoMin)
                  {
                     (*nitr)->clumpid = CLUMPONLIST;
                     clumpees.push_back(*nitr);
                  }
                  else
                     (*nitr)->clumpid = CLUMPNONE;

               // start a new clump.
               size_t nocp = 0;
               while (not clumpees.empty())
               {
                  partT* curf = clumpees.back();
                  clumpees.pop_back();

                  curf->clumpid = nclumpid;
                  nocp++;

                  neighs.clear();
                  neighs = NFW(curf, curf->h);

                  for (nitr = neighs.begin(); nitr != neighs.end(); nitr++)
                  {
                     if ((*nitr)->rho > _rhoMin)
                     {
                        if ((*nitr)->clumpid == nclumpid)
                           continue;

                        switch ((*nitr)->clumpid)
                        {
                        case CLUMPNONE:
                           (*nitr)->clumpid = CLUMPONLIST;
                           clumpees.push_back(*nitr);
                           break;

                        case CLUMPNOTSET:
                           (*nitr)->clumpid = CLUMPONLIST;
                           clumpees.push_back(*nitr);
                           break;

                        case CLUMPONLIST:
                           break;

                        default:
                           nclumpid = mergeClumps((*nitr)->clumpid,
                                                  nclumpid, _parts);
                           (*nitr)->clumpid = nclumpid;
                           break;
                        }
                     }
                     else
                     {
                        (*nitr)->clumpid = CLUMPNONE;
                     }
                  }
               }
            }
               break;

            case 1:
            {
               cPart->clumpid = *(neighClumps.begin());
            }
               break;

            default:
            {
               cType nclumpid = *(neighClumps.begin());

               std::set<cType>::const_iterator ncItr;
               for (ncItr = neighClumps.begin(); ncItr != neighClumps.end();
                    ncItr++)
                  nclumpid = mergeClumps(*ncItr, nclumpid, _parts);

               cPart->clumpid = nclumpid;
            }
               break;
            }
         }
         else
         { }
      }
   }

   Tree.clear();
}

template<typename _partT>
void Clumps<_partT>::potSearch(ParticleSet<_partT>& _parts)
{
   const size_t        nop           = _parts.getNop();
   const fType         G             = _parts.attributes["gravconst"];
   treeT::czllPtrVectT CZbottomLoc   = Tree.getCZbottomLoc();
   const int           noCZbottomLoc = CZbottomLoc.size();

   gravT gravWorker(&Tree, G);

#pragma omp parallel for firstprivate(gravWorker)
   for (int i = 0; i < noCZbottomLoc; i++)
      gravWorker.calcGravity(CZbottomLoc[i]);

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

template<typename _partT>
void Clumps<_partT>::potSearchOld(ParticleSet<_partT>& _parts)
{
   const size_t nop = _parts.getNop();
   const fType  G   = _parts.attributes["gravconst"];

   treeT::czllPtrVectT CZbottomLoc   = Tree.getCZbottomLoc();
   const int           noCZbottomLoc = CZbottomLoc.size();

   gravT gravWorker(&Tree, G);

#pragma omp parallel for firstprivate(gravWorker)
   for (int i = 0; i < noCZbottomLoc; i++)
      gravWorker.calcGravity(CZbottomLoc[i]);

   //NeighFindWorker<_partT> NFW(&Tree);
   typedef NeighFindSortedWorker<_partT>                     NFSWT;
   typedef typename NeighFindSortedWorker<_partT>::pptrDLT   partPtrDLT;
   NFSWT NFSW(&Tree);

   partPtrLT potRanked;

   for (size_t i = 0; i < nop; i++)
      potRanked.push_back(&(_parts[i]));

   lowerPot potSorter;
   potRanked.sort(potSorter);

   //typename partPtrLT::iterator pItr, nItr;
   typename partPtrLT::iterator pItr;
   typename partPtrDLT::iterator nItr;

   nextFreeId = 1;
   cType cID = nextFreeId;

   for (pItr = potRanked.begin(); pItr != potRanked.end(); pItr++)
   {
      partT       * cPart = *pItr;
      const fType cPot    = cPart->pot;
      const fType srad    = 2. * cPart->h;

      partPtrDLT neighs = NFSW(cPart, srad);


      if (neighs.size() == 0)
      {
         cPart->clumpid = CLUMPNONE;
         continue;
      }

      if (cPart->clumpid == CLUMPNOTSET)
      {
         cID            = nextFreeId;
         cPart->clumpid = cID;
         nextFreeId++;
      }
      else
         cID = cPart->clumpid;

      fType maxPot = -fTypeInf;
      for (nItr = neighs.begin(); nItr != neighs.end(); nItr++)
      {
         const fType nPot = (*nItr).ptr->pot;

         if (nPot < maxPot)
            break;

         if (nPot > cPot)
         {
            (*nItr).ptr->clumpid = cID;
            if (nPot > maxPot)
               maxPot = nPot;
         }
      }
   }
}

template<typename _partT>
void Clumps<_partT>::collectClumps(ParticleSet<_partT>& _parts,
                                   const fType          _minMass)
{
   const size_t nop = _parts.getNop();
   const size_t noc = nextFreeId;

   std::vector<clumpT> massRanked(noc);

   for (size_t i = 0; i < nop; i++)
   {
      if (_parts[i].clumpid < CLUMPNONE)
         _parts[i].clumpid = CLUMPNONE;

      const size_t cID = _parts[i].clumpid;
      assert(cID < noc);
      massRanked[cID].addParticle(&_parts[i]);
   }

   std::sort(massRanked.begin(), massRanked.end(), higherClumpMass());

   std::cout << "mass ranked     max: " << massRanked[0].m << "\n";

   std::vector<size_t> massRankedInv(noc);
   for (size_t i = 0; i < noc; i++)
      massRankedInv[massRanked[i].id] = i;


   size_t maxID = 0;
   for (size_t i = 0; i < noc; i++)
   {
      if (massRanked[i].m < _minMass)
      {
         std::cout << massRanked[i].m << " ";
         break;
      }
      maxID = i;
   }
   std::cout << "maxID " << maxID << "\n";

   for (size_t i = 0; i < nop; i++)
   {
      const size_t oldID = _parts[i].clumpid;
      const size_t newID = massRankedInv[oldID];

      if (newID <= maxID)
         _parts[i].clumpid = newID;
      else
         _parts[i].clumpid = CLUMPNONE;
   }

   const size_t fnoc = maxID + 1;

   parentT::resize(fnoc);

   for (size_t i = 0; i < nop; i++)
   {
      clumpT& curClump(parentT::operator[](_parts[i].clumpid));
      curClump.addParticle(&(_parts[i]));
   }
   std::cout << "added parts   \n";


   for (size_t i = 0; i < fnoc; i++)
   {
      clumpT& curClump(parentT::operator[](i));
      curClump.calcProperties();
      curClump.calcAngularMom();
      curClump.assignID(i + 1);
   }

   std::cout << "calc props    \n";

   return;
}
}

#endif
