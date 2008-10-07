#ifndef OOSPHEXECUTOR_H
#define OOSPHEXECUTOR_H

#include <boost/type_traits/remove_pointer.hpp>
#include <boost/type_traits/add_pointer.hpp>
#include <boost/mpl/for_each.hpp>

#include <boost/progress.hpp>

#include <iostream>

#include "manager.h"
#include "particle.h"
#include "costzone.h"

namespace oosph
{


namespace mpl = boost::mpl;

/**
 * \brief Driver Class for the Executor Hack
 * 
 * \author Pascal Bauer pbauer@phim.unibe.ch
 */


template <typename ExecutorType>
class Driver
{
public:

    template <typename CollectionType>
    inline void operator()(CollectionType) const
    {
        boost::remove_pointer<CollectionType>::type::Local(ExecutorType::Particle);
    }
protected:
private:
};


/**
 * \brief Prelocal Wraper for the Executor Hack.
 * 
 * \author Pascal Bauer pbauer@phim.unibe.ch
 */
class PreLocal
{
public:

    template <typename CollectionType>
    inline void operator()(CollectionType) const
    {
        boost::remove_pointer<CollectionType>::type::PreLocal();
    }
protected:
private:
};

/**
 * \brief PostLocal Wrapper for the Executor Hack.
 * 
 * \author Pascal Bauer pbauer@phim.unibe.ch
 */
class PostLocal
{
public:

    template <typename CollectionType>
    inline void operator()(CollectionType) const
    {
        boost::remove_pointer<CollectionType>::type::PostLocal();
    }
protected:
private:
};

/**
 * \brief Executor Class for the SPH-Particle Functions.
 * 
 * \author Pascal Bauer pbauer@phim.unibe.ch
 */

template <typename Collection, typename SPH_Trait>
class Executor
{
public:

    typedef Executor<Collection, SPH_Trait> self_type;
    typedef Executor<Collection, SPH_Trait>& self_reference;
    typedef Executor<Collection, SPH_Trait>* self_pointer;

    typedef SPH_Trait sph_trait;
    typedef typename sph_trait::simulation_trait sim_trait;

    typedef typename sph_trait::neighbour_type neighbour_type;
    typedef typename sph_trait::neighbour_reference neighbour_reference;
    typedef typename sph_trait::neighbour_pointer neighbour_pointer;

    typedef Manager<sim_trait>  manager_type;
    typedef Manager<sim_trait>& manager_reference;
    typedef Manager<sim_trait>* manager_pointer;

    typedef Driver<Executor<Collection, SPH_Trait> > driver_type;

    typedef mpl::vector_c<size_t, oosph::X> CostZoneIndex;
    typedef CostZone<CostZoneIndex, sim_trait> costzone_type;

    void operator()(void);
    static size_t Particle;

    static self_reference Instance(void);

    costzone_type& CostZone;

protected:

    Executor(void);
    ~Executor(void);

private:
    static manager_reference Manager;
    static neighbour_reference NB;
    static self_pointer _instance;

};

template <typename Collection, typename SPH_Trait>
typename Executor<Collection, SPH_Trait>::self_pointer
Executor<Collection, SPH_Trait>::_instance = NULL;

template <typename Collection, typename SPH_Trait>
typename Executor<Collection, SPH_Trait>::self_reference
Executor<Collection, SPH_Trait>::Instance(void)
{
    if (_instance != NULL)
        return *_instance;
    else
    {
        _instance = new self_type;
        return *_instance;
    }
}

template <typename Collection, typename SPH_Trait>
typename Executor<Collection, SPH_Trait>::manager_reference Executor<Collection, SPH_Trait>::Manager(manager_type::Instance());

template <typename Collection, typename SPH_Trait>
size_t Executor<Collection, SPH_Trait>::Particle = 0;

template <typename Collection, typename SPH_Trait>
typename Executor<Collection, SPH_Trait>::neighbour_reference
Executor<Collection, SPH_Trait>::NB(neighbour_type::Instance());

template <typename Collection, typename SPH_Trait>
Executor<Collection, SPH_Trait>::Executor(void) : CostZone(costzone_type::Instance())
{}
;

template <typename Collection, typename SPH_Trait>
Executor<Collection, SPH_Trait>::~Executor(void)
{}
;

template <typename Collection, typename SPH_Trait>
void Executor<Collection, SPH_Trait>::operator()(void)
{
    mpl::for_each<Collection>(PreLocal() );

    /**  \todo bring this in the typedef section */
    const size_t Size = SPH_Trait::manager_type::mem_manager_type::Instance().Data.size1();

    /** Start measuring time */
    CostZone.StartTimer();

    for (Particle = 0; Particle < Size; Particle++)
    {
        neighbour_type::GetNeighbours(Particle);
        mpl::for_each<Collection>(driver_type() );
    }

    /** Stop measuring time */
    CostZone.StopTimer();
    
    mpl::for_each<Collection>(PostLocal() );
};

}

#endif
