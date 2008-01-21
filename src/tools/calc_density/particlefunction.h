#ifndef OOSPHPARTICLEFUNCTION_H
#define OOSPHPARTICLEFUNCTION_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <cassert>

#include "manager.h"

namespace oosph
{

/**
 * \brief Parent Class for the SPH Particle Functions.
 * 
 * \author Pascal Bauer pbauer@phim.unibe.ch
 */

namespace num = boost::numeric::ublas;

template <typename T_LeafType, typename SPH_Trait>
class ParticleFunction
{
public:

    typedef ParticleFunction  self_type;
    typedef ParticleFunction& self_reference;
    typedef ParticleFunction* self_pointer;

    typedef T_LeafType  leaf_type;
    typedef T_LeafType& leaf_reference;
    typedef T_LeafType* leaf_pointer;

    typedef typename SPH_Trait::simulation_trait sim_trait;

    typedef typename SPH_Trait::neighbour_type      neighbour_type;
    typedef typename SPH_Trait::neighbour_reference neighbour_reference;
    typedef typename SPH_Trait::neighbour_pointer   neighbour_pointer;

    typedef typename sim_trait::matrix_type matrix_type;
    typedef typename sim_trait::matrix_reference matrix_reference;
    typedef typename sim_trait::matrix_pointer matrix_pointer;

    typedef typename neighbour_type::matrix_row matrix_row;

    typedef typename neighbour_type::particle_container_type      particle_container_type;
    typedef typename neighbour_type::particle_container_reference particle_container_reference;
    typedef typename neighbour_type::particle_container_pointer   particle_container_pointer;

    typedef typename particle_container_type::iterator particle_iterator;
    typedef typename particle_container_type::const_iterator const_particle_iterator;

    typedef Manager<sim_trait>  manager_type;
    typedef Manager<sim_trait>& manager_reference;
    typedef Manager<sim_trait>* manager_pointer;

    typedef typename manager_type::mem_manager_type mem_manager_type;
    typedef typename manager_type::mem_manager_type& mem_manager_reference;
    typedef typename manager_type::mem_manager_type* mem_manager_pointer;

    typedef typename SPH_Trait::simulation_trait::matrix_column matrix_column;
    typedef typename SPH_Trait::simulation_trait::matrix_row_range matrix_row_range;
    typedef typename SPH_Trait::simulation_trait::value_type scalar_type;

    typedef num::range range;

    ParticleFunction(void);
    ~ParticleFunction(void);

    leaf_reference AsLeaf(void);

    static void PreLocal(void);
    static void Local(const size_t& idx);
    static void PostLocal(void);

protected:

    static particle_container_reference NeighbourContainer;

private:
    static manager_reference Manager;

};

template <typename T_LeafType, typename SPH_Trait>
typename ParticleFunction<T_LeafType, SPH_Trait>::manager_reference ParticleFunction<T_LeafType, SPH_Trait>::Manager(manager_type::Instance() );

template <typename T_LeafType, typename SPH_Trait>
typename ParticleFunction<T_LeafType, SPH_Trait>::particle_container_reference ParticleFunction<T_LeafType, SPH_Trait>::NeighbourContainer(neighbour_type::neighbours);

template <typename T_LeafType, typename SPH_Trait>
ParticleFunction<T_LeafType, SPH_Trait>::ParticleFunction(void)
{}

template <typename T_LeafType, typename SPH_Trait>
ParticleFunction<T_LeafType, SPH_Trait>::~ParticleFunction(void)
{}

template <typename T_LeafType, typename SPH_Trait>
typename ParticleFunction<T_LeafType, SPH_Trait>::leaf_reference ParticleFunction<T_LeafType, SPH_Trait>::AsLeaf(void)
{
    return static_cast<leaf_reference>(*this);
}

template <typename T_LeafType, typename SPH_Trait>
inline void ParticleFunction<T_LeafType, SPH_Trait>::PreLocal(void)
{
    leaf_type::PreLocal();
}

template <typename T_LeafType, typename SPH_Trait>
inline void ParticleFunction<T_LeafType, SPH_Trait>::PostLocal(void)
{
    leaf_type::PostLocal();
}

template <typename T_LeafType, typename SPH_Trait>
inline void ParticleFunction<T_LeafType, SPH_Trait>::Local(const size_t& idx)
{

    static matrix_reference Data(SPH_Trait::manager_type::mem_manager_type::Instance().Data);

    assert(idx < Data.size1() );

    matrix_row Current(Data, idx);

    particle_iterator END = NeighbourContainer.end();
    
    for (particle_iterator Pos = NeighbourContainer.begin(); Pos != END; Pos++)
    {
        assert (Pos->index() < Pos->data().size1() );

        leaf_type::Sum_N(Current, *Pos);
    }
}


}

#endif
