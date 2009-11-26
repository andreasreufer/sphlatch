#ifndef SPHLATCH_SPH_ALGORITHMS_CPP
#define SPHLATCH_SPH_ALGORITHMS_CPP

#include "typedefs.h"

namespace sphlatch {
template<typename _partT>
struct densSum
{
   void zero(_partT* const _i)
   {
      std::cout << "zero!\n";

      _i->rho = 0.;
   }

   void operator()(const _partT* _i, const _partT* const _j)
   {
      const vect3dT rvec = _i->pos - _j->pos;
      const fType   rr   = rvec[0] * rvec[0] +
                           rvec[1] * rvec[1] +
                           rvec[2] * rvec[2];

      std::cout << sqrt(rr) << "\n";
      //_i->rho += 42.*(_j->m);
   }
};

/*template<typename _partT>
   void densSum<_partT>::zero(const _partT* _i)
   {
    _i->rho = 0.;
   };

   template<typename _partT>
   void densSum<_partT>::operator()(const _partT* _i)
   {
    const vect3dT rvec = _i->pos - _j->pos;
    std::cout << rvec << "\n";
   };*/
};

#endif
