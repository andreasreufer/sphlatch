#ifndef SPHLATCH_SPH_ALGORITHMS_CPP
#define SPHLATCH_SPH_ALGORITHMS_CPP

#include "typedefs.h"

namespace sphlatch {

template<typename _partT>
struct densSum
{
  void zero(const _partT* _i)
  {
    _i->rho = 0.;
  }

  void operator()(const _partT* _i, const _partT* const _j)
  {
    _i->rho += 42.*(_j->m);
  }
};

};

#endif
