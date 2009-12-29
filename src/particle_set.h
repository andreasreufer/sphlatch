#ifndef SPHLATCH_PARTICLE_SET_H
#define SPHLATCH_PARTICLE_SET_H

#include <vector>
#include <string>
#include <sys/stat.h>

#ifdef SPHLATCH_MPI
 #include <mpi.h>
#endif

#ifdef SPHLATCH_HDF5
 #define H5_NO_DEPRECATED_SYMBOLS
 #include <hdf5.h>
#endif

#include "typedefs.h"
#include "io_particle.h"

namespace sphlatch {
template<typename _partT>
class ParticleSet {
public:
   ParticleSet();
   ~ParticleSet();

   typedef std::list<std::string>   strLT;
   typedef partT::ioVar             ioVarT;
   typedef partT::ioVarLT           ioVarLT;
   typedef IOPart                   IOPartT;

   _partT& operator[](const size_t _i);
   void operator=(const ParticleSet& _ps);

#ifdef SPHLATCH_HDF5
   void saveHDF5(std::string _file);
   void loadHDF5(std::string _file);

   void singlePrecOut();
   void doublePrecOut();
#endif

   void resize(const size_t _i);
   size_t getNop();

   vect3dT getCom();
   box3dT  getBox();

   cType   step;

   ioVarLT loadVars, saveVars;
   attrMT  attributes;

protected:
   std::vector<_partT> parts;


#ifdef SPHLATCH_HDF5
   hid_t getLocFhRW(std::string _file);
   hid_t getLocFhRO(std::string _file);
   bool objExist(hid_t _fh, std::string _op);

   static herr_t getObsCB(hid_t _fh, const char* _name,
                          const H5L_info_t* info, void* opData);
   static herr_t getAttrsCB(hid_t _loc, const char* _name,
                            const H5A_info_t* info, void* opData);

   static strLT dsetsInFile, attrsInFile;
   hid_t        h5mITYPE, h5mFTYPE, h5fITYPE, h5fFTYPE;
#endif
};
};

#endif
