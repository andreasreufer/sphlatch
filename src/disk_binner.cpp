#ifndef SPHLATCH_DISK_BINNER_CPP
#define SPHLATCH_DISK_BINNER_CPP

/*
 *  disk_binner.cpp
 *
 *
 *  Created by Andreas Reufer on 17.10.10
 *  Copyright 2010 University of Berne. All rights reserved.
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

#include "fof.cpp"

namespace sphlatch {
class DiskBins {
public:
   DiskBins(const size_t _size)
   {
      Ltar.resize(_size, 3);
      Limp.resize(_size, 3);
      Ltot.resize(_size, 3);

      mtar.resize(_size);
      mimp.resize(_size);
      mtot.resize(_size);

      for (size_t i = 0; i < _size; i++)
      {
         mtar(i) = 0.;
         mimp(i) = 0.;
         mtot(i) = 0.;

         Ltar(i, 0) = 0.;
         Ltar(i, 1) = 0.;
         Ltar(i, 2) = 0.;

         Limp(i, 0) = 0.;
         Limp(i, 1) = 0.;
         Limp(i, 2) = 0.;

         Ltot(i, 0) = 0.;
         Ltot(i, 1) = 0.;
         Ltot(i, 2) = 0.;
      }
   }

   ~DiskBins() { }

   fvectT mtar, mimp, mtot;
   fmatrT Ltar, Limp, Ltot;

   void add(const size_t _i, fType _mi, vect3dT _Li, const bool _isTarg)
   {
      if (_isTarg)
      {
         mtar(_i) += _mi;

         Ltar(_i, 0) += _Li[0];
         Ltar(_i, 1) += _Li[1];
         Ltar(_i, 2) += _Li[2];
      }
      else
      {
         mimp(_i) += _mi;

         Limp(_i, 0) += _Li[0];
         Limp(_i, 1) += _Li[1];
         Limp(_i, 2) += _Li[2];
      }

      mtot(_i) += _mi;

      Ltot(_i, 0) += _Li[0];
      Ltot(_i, 1) += _Li[1];
      Ltot(_i, 2) += _Li[2];
   }

   void save(HDF5File& _hf, std::string _sfx)
   {
      _hf.savePrimitive("Ltar" + _sfx, Ltar);
      _hf.savePrimitive("Limp" + _sfx, Limp);
      _hf.savePrimitive("Ltot" + _sfx, Ltot);

      _hf.savePrimitive("mtar" + _sfx, mtar);
      _hf.savePrimitive("mimp" + _sfx, mimp);
      _hf.savePrimitive("mtot" + _sfx, mtot);
   }
};

template<typename _partT>
class DiskBinner {
public:
   typedef ParticleSet<_partT>   partSetT;

public:
   DiskBinner(partSetT& _parts) : parts(_parts), maxMat(32) { }
   ~DiskBinner() { }

   void saveBins(const iType _ccid, const fType _rmin, const fType _rmax,
                 const size_t _noBins, std::string _fname);
   void findFOF(const iType _ccid, const fType _rhoMin, const fType _hmult,
                const fType _mMin);

private:
   partSetT&    parts;
   const size_t maxMat;
};


template<typename _partT>
void DiskBinner<_partT>::saveBins(const iType _ccid, const fType _rmin,
                                  const fType _rmax, const size_t _noBins,
                                  std::string _fname)
{
   iType idfilt = 2000000;
   if (parts.attributes.count("idfilt") > 0)
      idfilt = parts.attributes["idfilt"];

   // determine center of mass of the clump particles
   // create new particle set, creating all particles NOT in the clump
   const size_t nop = parts.getNop();

   fType   mclmp = 0.;
   vect3dT rcom, vcom;
   rcom = 0., 0., 0.;
   vcom = 0., 0., 0.;

   for (size_t i = 0; i < nop; i++)
      if ((parts[i].clumpid == _ccid) && (parts[i].orbit == 1))
      {
         const fType m = parts[i].m;
         mclmp += m;
         rcom  += m * parts[i].pos;
         vcom  += m * parts[i].vel;
      }

   rcom /= mclmp;
   vcom /= mclmp;

   vect3dT Lcom;
   Lcom = 0., 0., 0.;

   for (size_t i = 0; i < nop; i++)
      if ((parts[i].clumpid == _ccid) && (parts[i].orbit == 1))
      {
         const vect3dT rrel = rcom - parts[i].pos;
         const vect3dT vrel = vcom - parts[i].vel;

         Lcom += cross(rrel, vrel) * parts[i].m;
      }

   std::vector<DiskBins*> matbins;
   DiskBins totbins(_noBins);

   matbins.resize(maxMat);
   for (size_t i = 0; i < maxMat; i++)
      matbins[i] = NULL;

   fvectT r;
   r.resize(_noBins);

   const fType dr = (_rmax - _rmin) / (_noBins);
   for (size_t i = 0; i < _noBins; i++)
      r(i) = (0.5 + static_cast<fType>(i)) * dr;

   // make sure this materials have disk bins, even if
   // no such particles are present. this makes plotting easier
   matbins[0] = new DiskBins(_noBins);
   matbins[1] = new DiskBins(_noBins);
   matbins[2] = new DiskBins(_noBins);
   matbins[4] = new DiskBins(_noBins);
   matbins[5] = new DiskBins(_noBins);

   for (size_t i = 0; i < nop; i++)
      if (parts[i].clumpid == _ccid)
      {
         const iType mati = parts[i].mat;

         if (matbins[mati] == NULL)
            matbins[mati] = new DiskBins(_noBins);

         const vect3dT rrel = rcom - parts[i].pos;
         const vect3dT vrel = vcom - parts[i].vel;

         const fType   mi  = parts[i].m;
         const fType   ri  = sqrt(dot(rrel, rrel));
         const vect3dT Li  = cross(rrel, vrel) * parts[i].m;
         const size_t  idx = std::min(static_cast<long int>(_noBins - 1),
                                      lrint(ri / dr));

         matbins[mati]->add(idx, mi, Li, (parts[i].id < idfilt));
         totbins.add(idx, mi, Li, (parts[i].id < idfilt));
      }

   std::string nostr   = boost::lexical_cast<std::string>(_ccid);
   std::string diskstr = "disk";
   for (size_t i = nostr.size(); i < 3; i++)
      diskstr.append("0");
   diskstr.append(nostr);

   HDF5File diskFile(_fname);
   diskFile.doublePrecOut();

   // save attributes of particle dump to disk data
   diskFile.setNewRoot(parts.getStepName());
   diskFile.saveAttributes(parts.attributes);

   // now save the disk data
   diskFile.setNewRoot(parts.getStepName() + "/" + diskstr);
   diskFile.savePrimitive("r", r);

   for (size_t i = 0; i < maxMat; i++)
      if (matbins[i] != NULL)
      {
         std::stringstream matsfx, idstr;
         idstr << i;
         matsfx << "_mat";
         for (size_t j = idstr.str().size(); j < 2; j++)
            matsfx << "0";
         matsfx << i;
         matbins[i]->save(diskFile, matsfx.str());
      }
   totbins.save(diskFile, "_tot");
}

template<typename _partT>
void DiskBinner<_partT>::findFOF(const iType _ccid, const fType _rhoMin,
                                 const fType _hmult,
                                 const fType _mMin)
{
   const size_t notp = parts.getNop();
   partSetT     disk;

   for (size_t i = 0; i < notp; i++)
      if (parts[i].clumpid == _ccid)
         (disk.insert(parts[i])).id = i;

   const size_t nodp = disk.getNop();

   // search for individual clumps
   // separate non-singleton Tree
   treeT Tree;
   Tree.setExtent(disk.getBox() * 1.1);

   for (size_t i = 0; i < nodp; i++)
      Tree.insertPart(disk[i]);

   Tree.update(0.8, 1.2);

   for (size_t i = 0; i < nodp; i++)
      disk[i].friendid = sphlatch::FRIENDNOTSET;

   sphlatch::FOF<partT> fof(disk, Tree);
   fof.search(_rhoMin, _hmult, _mMin);

   for (size_t i = 0; i < nodp; i++)
      parts[disk[i].id].friendid = disk[i].friendid;

   Tree.clear();
}
}
#endif
