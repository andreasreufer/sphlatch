#ifndef SPH_IO_PARTICLE_H
#define SPH_IO_PARTICLE_H

/*
 *  sph_io_particle.h
 *
 *
 *  Created by Andreas Reufer on 12.11.09
 *  Copyright 2009 University of Berne. All rights reserved.
 *
 */

#include <stddef.h>
#include "typedefs.h"


namespace sphlatch {
class IOPart {
public:
   enum storageType { FTYPE, ITYPE };

   class ioVar {
public:
      ioVar(const size_t _off, const size_t _wdt, storageType _type)
      {
         offset = _off;
         width  = _wdt;
         type   = _type;
      }

      size_t      offset, width;
      storageType type;
   };

   typedef std::list<ioVar>   ioVarLT;


   ioVar storeVar(const fType& _f, std::string _name)
   {
      const size_t offset = reinterpret_cast<const char*>(&_f) -
                            reinterpret_cast<char*>(this);

      return(ioVar(offset, 1, FTYPE));
   }

   ioVar storeVar(const iType& _i, std::string _name)
   {
      const size_t offset = reinterpret_cast<const char*>(&_i) -
                            reinterpret_cast<char*>(this);

      return(ioVar(offset, 1, ITYPE));
   }

   ioVar storeVar(const vect3dT& _v, std::string _name)
   {
      const size_t offset = reinterpret_cast<const char*>(&_v[0]) -
                            reinterpret_cast<char*>(this);

      return(ioVar(offset, 3, FTYPE));
   }
};
};
#endif
