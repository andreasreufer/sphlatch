#ifndef SPHLATCH_PARTICLE_SET_H
#define SPHLATCH_PARTICLE_SET_H

#include <vector>
#include <string>
#include <sys/stat.h>

#ifdef SPHLATCH_MPI
 #include <mpi.h>
#endif

#ifdef SPHLATCH_HDF5
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

   _partT& operator[](const size_t _i);

#ifdef SPHLATCH_HDF5
   void saveDump(std::string _file);
   void loadDump(std::string _file);
#endif

   void resize(const size_t _i);

   cType step;
   typedef std::list<std::string> strLT;
   typedef partT::ioVar ioVarT;
   typedef partT::ioVarLT ioVarLT;
   typedef IOPart IOPartT;

   ioVarLT loadVars, saveVars;

private:
   std::vector<_partT> parts;


#ifdef SPHLATCH_HDF5
   hid_t getLocFhRW(std::string _file);
   hid_t getLocFhRO(std::string _file);

   static herr_t getObsCB(hid_t _fh, const char* _name,
                          const H5L_info_t* info, void* opData);
   static herr_t getAttrsCB(hid_t _loc, const char* _name,
                            const H5A_info_t* info, void* opData);

   static strLT dsetsInFile, attrsInFile;
#endif
};
};

#endif
