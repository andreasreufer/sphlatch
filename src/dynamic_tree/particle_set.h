#ifndef SPHLATCH_PARTICLE_SET_H
#define SPHLATCH_PARTICLE_SET_H

#include <vector>
#include <string>
#include <sys/stat.h>

#ifdef SPHLATCH_MPI
#include <mpi.h>
#endif

#include <hdf5.h>

#include "typedefs.h"

namespace sphlatch {

template<typename _partT>
class ParticleSet {
  public:

  _partT&  operator[](const size_t _i);

  void saveDump(std::string _file);
  void loadDump(std::string _file);

  void resize(const size_t _i);


  private:
      std::vector<_partT> parts;

      /*hid_t getLocFhRW(std::string _file);
      hid_t getLocFhRO(std::string _file);
      static herr_t getObsCB();*/

};
};

#endif
