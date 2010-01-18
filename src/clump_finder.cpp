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

#include "typedefs.h"
#include "particle_set.cpp"
#include "clump_particle.h"

#include "bhtree.cpp"
#include "bhtree_worker_neighfunc.cpp"

namespace sphlatch {
template<typename _partT>
class Clumps
{
public:
   typedef _partT                partT;
   typedef ParticleSet<_partT>   partSetT;
   typedef std::list<_partT*>    partPtrLT;

   typedef BHTree                treeT;

public:
   Clumps() { }
   ~Clumps() { }

   void findClumps(partSetT& _parts, const fType _rhoMin);
   void startNewClump(partT* _part);
   cType mergeClumps(const cType cida, const cType cidb, partSetT& _parts);
};

template<typename _partT>
cType Clumps<_partT>::mergeClumps(const cType cida, const cType cidb,
                                  partSetT& _parts)
{
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
void Clumps<_partT>::findClumps(ParticleSet<_partT>& _parts,
                                const fType          _rhoMin)
{
   treeT& Tree(treeT::instance());

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

   cType nextFreeId = nop + 1;

   NeighFindWorker<_partT> NFW(&Tree);
   for (size_t i = 0; i < nop; i++)
   //for (size_t i = 0; i < 1; i++)
   {
      //std::cout << i << "\n";
      //std::cout << "particle " << i << " rho = " << _parts[i].rho << "\n";
      cType& clumpid(_parts[i].clumpid);
      partT  * cPart = &(_parts[i]);

      //
      // if density is too low, particle does not belong to any clump
      // regardless of previous state
      //
      if (_parts[i].rho < _rhoMin)
      {
         clumpid = CLUMPNONE;
         //std::cout << " rho too low, continue ...\n";
         continue;
      }
      else
      {
         if ((clumpid == CLUMPNONE) || (clumpid == CLUMPNOTSET))
         {
            partPtrLT friends = NFW(cPart, cPart->h);

            std::set<cType> neighClumps;
            typename partPtrLT::iterator fitr = friends.begin();
            while (fitr != friends.end())
            {
               if ((*fitr)->clumpid > 0)
                  neighClumps.insert((*fitr)->clumpid);
               fitr++;
            }

            const size_t noClumpsFound = neighClumps.size();
            switch (noClumpsFound)
            {
            case 0:
            {
               clumpid = nextFreeId;
               nextFreeId++;

               //std::cout << "newclump! " << clumpid << " " << friends.size() << "\n";

               while (friends.size() > 0)
               {
                  partT* curf = friends.back();
                  friends.pop_back();

                  if (curf->rho > _rhoMin)
                  {
                     if (curf->clumpid < 1)
                     {
                     //curf->clumpid = clumpid;
                        //std::cout << " untouched!\n";
                     }
                     else
                     {
                        //FIXME: merge clumps
                        //std::cout << " was " << curf->clumpid << "\n";
                        //clumpid = mergeClumps(curf->clumpid, clumpid, _parts);
                     }
                     curf->clumpid = clumpid;
                     //std::cout << " flavoured!\n";
                  }
               }
               //startNewClump(cPart);
            }
               break;

            case 1:
            {
               cPart->clumpid = *(neighClumps.begin());
            }
               break;

            default:
            {
               //FIXME: merge clumps
               cPart->clumpid = *(neighClumps.begin());
            }
               break;
            }
         }
         else
         { }
      }

      //
      // density is high enough to belong to a clump
      //

      /*if ((clumpid == CLUMPNONE) || (clumpid == CLUMPNOTSET))
         {

         std::cout << " rho high enough, new clump! (id: "
                   << nextFreeId << ")\n";
         clumpid = nextFreeId;
         nextFreeId++;
         }*/

      //
      // now start to look for friends
      //

      /*std::list<_partT*> ffriends;
         ffriends.push_back(&(_parts[i]));

         while (ffriends.size() != 0)
         {
         //std::cout << " ff: " << ffriends.size() << "\n";
         _partT* const curff = ffriends.back();
         ffriends.pop_back();

         if ((ffriends.size() % 1000) == 0)
            std::cout << ffriends.size() << "\n";

         if (curff->rho > _rhoMin)
         {
            curff->clumpid = clumpid;
            // if two clumps are connected, create a single
            // new one with the lower id of both
            if ((curff->clumpid > 0) && (curff->clumpid != clumpid))
            {
               cType oldid, newid;
               if (clumpid > curff->clumpid)
               {
                  oldid = clumpid;
                  newid = curff->clumpid;
               }
               else
               {
                  oldid = curff->clumpid;
                  newid = clumpid;
               }

               std::cout << "clumps " << newid
                         << " & " << oldid
                         << " -> " << newid << "\n";

               for (size_t i = 0; i < nop; i++)
                  if (_parts[i].clumpid == oldid)
                     _parts[i].clumpid = newid;

               clumpid = newid;
            }


            const fType srad = std::min(curff->h, 3.e8);

            // add neighbours to the friends of friends list
            //std::list<_partT*> newffriends = NFW(curff, curff->h);
            std::list<_partT*> newffriends = NFW(curff, srad);
            //std::cout << "   " << newffriends.size() << " new friends found\n";

            typename std::list<_partT*>::iterator fitr = newffriends.begin();
            while (fitr != newffriends.end())
            {
               // only add it to the search, if
               if (((*fitr)->clumpid == -1) && ((*fitr)->rho > _rhoMin))
               {
                  ffriends.push_back(*fitr);
               }

               fitr++;
            }

            // only add the particle, if it does not yet belong to a clump
            //ffriends.splice(ffriends.end(), newffriends);
         }

         // condition for a friend
         //if curff->rho > rhoMin;
         }*/
   }



   //std::set<cType, cType> clumpSizes;
   std::set<cType> clumpIDs;

   // assign real clump IDs
   for (size_t i = 0; i < nop; i++)
   {
      clumpIDs.insert(_parts[i].clumpid);

      //if (_parts[i].clumpid > 0)
      //  std::cout << i << ": " << _parts[i].clumpid << "\n";
   }

   /*std::set<cType>::const_iterator ciditr;
      for (ciditr = clumpIDs.begin(); ciditr != clumpIDs.end(); ciditr++)
      std::cout << *ciditr << "\n";*/

   std::cout << clumpIDs.size() << " clumps known!\n";


   // output clumps?

   // clear tree
   Tree.clear();
}

template<typename _partT>
void Clumps<_partT>::startNewClump(partT* _part)
{ }
}

#endif
