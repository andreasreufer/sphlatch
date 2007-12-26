#ifndef OOSPH_NEIGHBOURS_H
#define OOSPH_NEIGHBOURS_H

#include <vector>
#include "manager.h"

namespace oosph
{

/**
 * \brief Parentclass for all the Neighbour Search Algorithms.
 * 
 * \author Pascal Bauer
 * 
 * This is the Parent class that should be used for any Neighbour Search Algorithm. It defines a minimal Interface, that will be used by the ParticleFunctiona Class
 */

template <typename T_LeafType, typename SimTrait>
class Neighbours
{
public:

    template <typename T1, typename T2>
    friend class ParticleFunction;

    typedef Neighbours  self_type;
    typedef Neighbours& self_reference;
    typedef Neighbours* self_pointer;

    typedef T_LeafType  leaf_type;
    typedef T_LeafType& leaf_reference;
    typedef T_LeafType* leaf_pointer;

    typedef typename SimTrait::matrix_row  matrix_row;

    typedef typename SimTrait::particle_container_type  particle_container_type;
    typedef typename SimTrait::particle_container_reference particle_container_reference;
    typedef typename SimTrait::particle_container_pointer particle_container_pointer;

    typedef Manager<SimTrait>  manager_type;
    typedef Manager<SimTrait>& manager_reference;
    typedef Manager<SimTrait>* manager_pointer;

    leaf_reference AsLeaf(void);

    static leaf_reference Construct(void);
    static particle_container_reference GetNeighbours(const size_t& index);


protected:

    static particle_container_reference neighbours;

private:

};

template <typename T_LeafType, typename SimTrait>
inline typename Neighbours<T_LeafType, SimTrait>::leaf_reference Neighbours<T_LeafType, SimTrait>::AsLeaf(void)
{
    return static_cast<leaf_reference>(*this);
}

template <typename T_LeafType, typename SimTrait>
inline typename Neighbours<T_LeafType, SimTrait>::leaf_reference Neighbours<T_LeafType, SimTrait>::Construct(void)
{
    leaf_type::Construct();
    return AsLeaf();
}

template <typename T_LeafType, typename SimTrait>
inline typename Neighbours<T_LeafType, SimTrait>::particle_container_reference Neighbours<T_LeafType, SimTrait>::GetNeighbours(const size_t& index)
{
    leaf_type::operator[](index);
}

template <typename T_LeafType, typename SimTrait>
typename Neighbours<T_LeafType, SimTrait>::particle_container_reference Neighbours<T_LeafType, SimTrait>::neighbours(leaf_type::neighbours);

}

#endif
