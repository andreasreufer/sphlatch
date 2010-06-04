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
#include "io_particle.h"

#include "bhtree.cpp"
#include "bhtree_worker_neighfunc.cpp"
#include "bhtree_worker_grav.cpp"

namespace sphlatch {
template<typename _partT>
class Clump : public IOPart
{
public:
   typedef std::list<_partT*>   partPtrLT;
   vect3dT pos;
   vect3dT vel;
   vect3dT L;
   vect3dT P;
   fType   m, rc, V;

   cType id, nop;

   void addParticle(_partT* _partPtr)
   {
      pList.push_back(_partPtr);
   }

   bool operator<(const Clump& _rhs)
   {
      return(m < _rhs.m);
   }

   void calcProperties()
   {
      typename partPtrLT::const_iterator pItr;
      m   = 0.;
      pos = 0., 0., 0.;
      V   = 0.;

      for (pItr = pList.begin(); pItr != pList.end(); pItr++)
      {
         const fType mi = (*pItr)->m;
         m   += mi;
         P   += mi * (*pItr)->vel;
         pos += mi * (*pItr)->pos;
         V   += mi / (*pItr)->rho;
      }
      pos /= m;
      vel  = P / m;
      rc   = pow(0.75 * V / M_PI, 1. / 3.);

      L   = 0., 0., 0.;
      nop = 0;
      for (pItr = pList.begin(); pItr != pList.end(); pItr++)
      {
         const fType   mi   = (*pItr)->m;
         const vect3dT rvec = (*pItr)->pos - pos;
         L += mi * cross(rvec, (*pItr)->vel);
         nop++;
      }
   }

   ioVarLT getSaveVars()
   {
      ioVarLT vars;

      vars.push_back(storeVar(pos, "pos"));
      vars.push_back(storeVar(vel, "vel"));
      vars.push_back(storeVar(L, "L"));
      vars.push_back(storeVar(P, "P"));

      vars.push_back(storeVar(m, "m"));
      vars.push_back(storeVar(rc, "rc"));

      vars.push_back(storeVar(id, "id"));
      vars.push_back(storeVar(nop, "nop"));

      return(vars);
   }

private:
   partPtrLT pList;
};


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

   class lowerPot
   {
public:
      bool operator()(_partT* _p1, _partT* _p2)
      {
         return(_p1->pot < _p2->pot);
      }
   };

public:
   class SimpleClump {
public:
      /*bool operator()(SimpleClump& _c1, SimpleClump& _c2)
      {
         return(_c1.m < _c2.m);
      }*/
   
      bool operator<(SimpleClump& _rhs)
      {
         return(m < _rhs.m);
      }

      cType id;
      fType m;
   };

   /*class SimpleClumpSorter {
   };*/


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
   //finalize(_parts);
   //

   const fType minMass = 1.1926e+24;
   collectClumps(_parts, minMass);
}

template<typename _partT>
void Clumps<_partT>::getClumpsFOF(ParticleSet<_partT>& _parts,
                                  const fType          _rhoMin,
                                  const fType          _hMult)
{
   prepareSearch(_parts);
   FOF(_parts, _rhoMin, _hMult);
   //finalize(_parts);
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


   std::cout << potRanked.size() << "\n";
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

   std::vector<SimpleClump> massRanked(noc);

   for (size_t i = 0; i < noc; i++)
   {
      massRanked[i].m  = 0.;
      massRanked[i].id = i;
   }

   for (size_t i = 0; i < nop; i++)
   {
      const size_t cID = _parts[i].clumpid;
      assert(cID < noc);

      massRanked[cID].m += _parts[i].m;
   }

   std::sort(massRanked.begin(), massRanked.end());
   //massRanked.sort();

   std::cout << "mass ranked     max: " << massRanked[0].m << "\n";

   std::vector<size_t> oldID2newID(noc);
   for (size_t i = 0; i < noc; i++)
   {
      oldID2newID[massRanked[i].id] = i;
   }

   for (size_t i = 0; i < nop; i++)
   {
      const size_t oldID = _parts[i].clumpid;
      _parts[i].clumpid = oldID2newID[oldID];
   }

   size_t maxId = 0;
   for (size_t i = 0; i < noc; i++)
   {
      if (massRanked[i].m < _minMass)
         break;
      maxId = i;
   }

   std::cout << "maxId: " << maxId << "\n";

   const size_t fnoc = maxId + 1;

   parentT::resize(fnoc);

   for (size_t i = 0; i < nop; i++)
   {
      if (static_cast<size_t>(_parts[i].clumpid) > maxId)
         _parts[i].clumpid = CLUMPNONE;
      else
      {
         clumpT& curClump(parentT::operator[](i));
         curClump.addParticle(&(_parts[i]));
      }
   }
   
   for (size_t i = 0; i < noc; i++)
   {
      clumpT& curClump(parentT::operator[](i));
      curClump.calcProperties();
   }
}

/*template<typename _partT>
   void Clumps<_partT>::calcClumps(ParticleSet<_partT>& _parts)
   {
   const size_t nop = _parts.getNop();

   parentT::resize(noc);

   std::set<cType>::const_iterator cidItr = clumpIDs.begin();

   for (size_t i = 0; i < noc; i++)
   {
      clumpT&     curClump(parentT::operator[](i));
      const cType cid = *cidItr;
      curClump.id = cid;

      for (size_t j = 0; j < nop; j++)
         if (_parts[j].clumpid == cid)
            curClump.addParticle(&(_parts[j]));

      curClump.calcProperties();

      cidItr++;
   }
   return;
   }*/
}

#endif
