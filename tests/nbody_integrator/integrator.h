#ifndef OOSPHINTEGRATOR_H
#define OOSPHINTEGRATOR_H

#include <cstdlib>

#include "manager.h"

namespace oosph
{

/**
\brief Integrator
*/

template <typename T_LeafType, typename SimTrait>
class Integrator
{
public:
    typedef T_LeafType leaf_type;
    typedef T_LeafType& leaf_reference;
    typedef T_LeafType* leaf_pointer;

    leaf_reference AsLeaf(void)
    {
        return static_cast<leaf_reference>(*this);
    }

    typedef typename SimTrait::value_type value_type;

    static value_type CFL(void);
    static void Predictor(value_type dt);
    static void Corrector(value_type dt);
    static void Step(value_type dt);

    typedef Manager<SimTrait>  manager_type;
    typedef Manager<SimTrait>& manager_reference;
    typedef Manager<SimTrait>* manager_pointer;

protected:

private:

};

template <typename T_LeafType, typename SimTrait>
typename SimTrait::value_type
Integrator<T_LeafType, SimTrait>::CFL(void)
{
    return leaf_type::CLF();
}

template <typename T_LeafType, typename SimTrait>
void Integrator<T_LeafType, SimTrait>::Predictor(value_type dt)
{
    return leaf_type::Predictor(dt);
}

template <typename T_LeafType, typename SimTrait>
void Integrator<T_LeafType, SimTrait>::Corrector(value_type dt)
{
    return leaf_type::Corrector(dt);
}

template <typename T_LeafType, typename SimTrait>
void Integrator<T_LeafType, SimTrait>::Step(value_type dt)
{
    return  leaf_type::Step(dt);
}

// **********************************************************************************************************



}

#endif
