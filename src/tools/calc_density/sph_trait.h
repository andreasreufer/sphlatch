#ifndef OOSPHSPH_TRAIT_H
#define OOSPHSPH_TRAIT_H

#include "manager.h"

#include "simulation_trait.h"

namespace oosph
{

/**
 * \brief Collection of Classes for the SPH-Calculations.
 * 
 * \author Pascal Bauer pbauer@phim.unibe.ch
 * 
 */

template <typename Kernel_Type, typename Viscosity_Type, typename EOS_Type, typename Neighbour_Type, typename Simulation_Trait = oosph::SimulationTrait<> >
class SPH_Trait
{
public:

    typedef Simulation_Trait simulation_trait;
    
    typedef Manager<simulation_trait>  manager_type;
    typedef Manager<simulation_trait>& manager_reference;
    typedef Manager<simulation_trait>* manager_pointer;

    typedef Kernel_Type  kernel_type;
    typedef Kernel_Type& kernel_reference;
    typedef Kernel_Type* kernel_pointer;

    typedef Viscosity_Type  viscosity_type;
    typedef Viscosity_Type& viscosity_reference;
    typedef Viscosity_Type* viscosity_pointer;

    typedef EOS_Type  eos_type;
    typedef EOS_Type& eos_reference;
    typedef EOS_Type* eos_pointer;

    typedef Neighbour_Type  neighbour_type;
    typedef Neighbour_Type& neighbour_reference;
    typedef Neighbour_Type* neighbour_pointer;

};

}

#endif
