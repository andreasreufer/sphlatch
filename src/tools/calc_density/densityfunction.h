#ifndef OOSPHDENSITYFUNCTION_H
#define OOSPHDENSITYFUNCTION_H

#include <cstdlib>
#include <cassert>

#include <boost/mpl/at.hpp>
#include <boost/mpl/vector_c.hpp>

#include "particlefunction.h"

namespace oosph
{

/**
 * 
 * \brief DensityFunction
 * 
 * \author Pascal Bauer pbauer@phim.unibe.ch
 */
template <typename Indices, typename SPH_Trait>
class DensityFunction : public ParticleFunction<DensityFunction<Indices, SPH_Trait>, SPH_Trait>
{
public:

    enum Index { M = mpl::at_c<Indices, 0>::type::value,
                 RHO = mpl::at_c<Indices, 1>::type::value};

    typedef DensityFunction  self_type;
    typedef DensityFunction& self_reference;
    typedef DensityFunction* self_pointer;

    typedef ParticleFunction<self_type, SPH_Trait>  parent_type;
    typedef ParticleFunction<self_type, SPH_Trait>& parent_reference;
    typedef ParticleFunction<self_type, SPH_Trait>* parent_pointer;

    typedef typename parent_type::manager_type manager_type;
    typedef typename parent_type::manager_reference manager_reference;

    typedef SPH_Trait sph_trait;

    typedef typename sph_trait::kernel_type      kernel_type;
    typedef typename sph_trait::kernel_reference kernel_reference;

    typedef typename sph_trait::simulation_trait::matrix_column matrix_column;
    typedef typename sph_trait::simulation_trait::zero_vector zero_vector;
    typedef typename sph_trait::neighbour_type::matrix_row matrix_row;

    typedef typename sph_trait::simulation_trait sim_trait;

    typedef typename sim_trait::matrix_type matrix_type;
    typedef typename sim_trait::matrix_reference matrix_reference;
    typedef typename sim_trait::matrix_pointer matrix_pointer;

    static self_reference Instance(void);

    static void PreLocal(void);
    static void Sum_N(matrix_row& Current, const matrix_row& Neighbour);
    static void PostLocal(void);

protected:

    DensityFunction(void);
    ~DensityFunction(void);

private:

    static self_pointer _instance;
    static kernel_reference Kernel;
    static manager_reference Manager;
};

template <typename Indices, typename SPH_Trait>
typename DensityFunction<Indices, SPH_Trait>::manager_reference DensityFunction<Indices, SPH_Trait>::Manager(manager_type::Instance() );

template <typename Indices, typename SPH_Trait>
typename DensityFunction<Indices, SPH_Trait>::kernel_reference DensityFunction<Indices, SPH_Trait>::Kernel(kernel_type::Instance() );

template <typename Indices, typename SPH_Trait>
typename DensityFunction<Indices, SPH_Trait>::self_pointer DensityFunction<Indices, SPH_Trait>::_instance = NULL;

template <typename Indices, typename SPH_Trait>
typename DensityFunction<Indices, SPH_Trait>::self_reference DensityFunction<Indices, SPH_Trait>::Instance(void)
{
    if (_instance != NULL)
        return *_instance;
    else
    {
        _instance = new DensityFunction;
        return *_instance;
    }
};

template <typename Indices, typename SPH_Trait>
DensityFunction<Indices, SPH_Trait>::DensityFunction(void)
{}

template <typename Indices, typename SPH_Trait>
DensityFunction<Indices, SPH_Trait>::~DensityFunction(void)
{}

template <typename Indices, typename SPH_Trait>
inline void DensityFunction<Indices, SPH_Trait>::PreLocal(void)
{}

template <typename Indices, typename SPH_Trait>
inline void DensityFunction<Indices, SPH_Trait>::PostLocal(void)
{
#ifndef NDEBUG
    matrix_reference Data(SPH_Trait::manager_type::mem_manager_type::Instance().Data);

    for (size_t i = 0; i < Data.size1(); ++i )
    {
        assert( Data(i, RHO) > 0 );
    }
#endif
}

template <typename Indices, typename SPH_Trait>
inline void DensityFunction<Indices, SPH_Trait>::Sum_N(matrix_row& Current, const matrix_row& Neighbour)
{
    assert( Kernel.K(Current, Neighbour) >= 0 );
    assert( Neighbour(M) >= 0. );

    Current(RHO) += Neighbour(M) * Kernel.K(Current, Neighbour);
}

}

#endif
