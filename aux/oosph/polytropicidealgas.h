#ifndef OOSPHPOLYTROPICIDEALGAS_H
#define OOSPHPOLYTROPICIDEALGAS_H

#include <cstdlib>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/mpl/at.hpp>

#include "eos.h"
#include "manager.h"

namespace oosph
{

namespace mpl = boost::mpl;
namespace num = boost::numeric::ublas;
/**
 * \brief Polytropic Ideal Gas
 * 
 * \author Pascal Bauer pbauer@phim.unibe.ch
 */

template <typename Indices, typename Sim_Trait>
class PolytropicIdealGas : public EOS<PolytropicIdealGas<Indices, Sim_Trait>, Sim_Trait>
{
public:

    enum Index { RHO = mpl::at_c<Indices, 0>::type::value,
                 E   = mpl::at_c<Indices, 1>::type::value,
                 P   = mpl::at_c<Indices, 2>::type::value};


    typedef PolytropicIdealGas  self_type;
    typedef PolytropicIdealGas& self_reference;
    typedef PolytropicIdealGas* self_pointer;

    typedef EOS<PolytropicIdealGas<Indices, Sim_Trait>, Sim_Trait> parent_type;
    typedef EOS<PolytropicIdealGas<Indices, Sim_Trait>, Sim_Trait>& parent_reference;
    typedef EOS<PolytropicIdealGas<Indices, Sim_Trait>, Sim_Trait>* parent_pointer;

    typedef Manager<Sim_Trait>  manager_type;
    typedef Manager<Sim_Trait>& manager_reference;
    typedef Manager<Sim_Trait>* manager_pointer;

    typedef typename manager_type::mem_manager_type mem_manager_type;

    typedef typename Sim_Trait::matrix_type matrix_type;
    typedef typename Sim_Trait::matrix_reference matrix_reference;

    typedef typename parent_type::matrix_column matrix_column;
    typedef typename parent_type::value_type value_type;

    static void Pressure(void);

    static self_reference Instance(void);

protected:

    PolytropicIdealGas(void);
    ~PolytropicIdealGas(void);

private:

    static matrix_reference Data;
    static matrix_reference GData;

    static self_pointer _instance;

    static manager_reference Manager;

    static const matrix_column rho;
    static const matrix_column grho;

    static const matrix_column e;
    static const matrix_column ge;

    static matrix_column p;
    static matrix_column gp;

    static value_type gamma;
};

template <typename Indices, typename Sim_Trait>
typename PolytropicIdealGas<Indices, Sim_Trait>::matrix_reference
PolytropicIdealGas<Indices, Sim_Trait>::Data(mem_manager_type::Instance().Data);

template <typename Indices, typename Sim_Trait>
typename PolytropicIdealGas<Indices, Sim_Trait>::matrix_reference
PolytropicIdealGas<Indices, Sim_Trait>::GData(mem_manager_type::Instance().GData);

template <typename Indices, typename Sim_Trait>
typename PolytropicIdealGas<Indices, Sim_Trait>::value_type PolytropicIdealGas<Indices, Sim_Trait>::gamma;

template <typename Indices, typename Sim_Trait>
const typename PolytropicIdealGas<Indices, Sim_Trait>::matrix_column PolytropicIdealGas<Indices, Sim_Trait>::rho(Data, RHO);

template <typename Indices, typename Sim_Trait>
const typename PolytropicIdealGas<Indices, Sim_Trait>::matrix_column PolytropicIdealGas<Indices, Sim_Trait>::grho(GData, RHO);

template <typename Indices, typename Sim_Trait>
const typename PolytropicIdealGas<Indices, Sim_Trait>::matrix_column PolytropicIdealGas<Indices, Sim_Trait>::e(Data, E);

template <typename Indices, typename Sim_Trait>
const typename PolytropicIdealGas<Indices, Sim_Trait>::matrix_column PolytropicIdealGas<Indices, Sim_Trait>::ge(GData, E);

template <typename Indices, typename Sim_Trait>
typename PolytropicIdealGas<Indices, Sim_Trait>::matrix_column PolytropicIdealGas<Indices, Sim_Trait>::p(Data, P);

template <typename Indices, typename Sim_Trait>
typename PolytropicIdealGas<Indices, Sim_Trait>::matrix_column PolytropicIdealGas<Indices, Sim_Trait>::gp(GData, P);

template <typename Indices, typename Sim_Trait>
typename PolytropicIdealGas<Indices, Sim_Trait>::manager_reference PolytropicIdealGas<Indices, Sim_Trait>::Manager(manager_type::Instance() );

template <typename Indices, typename Sim_Trait>
typename PolytropicIdealGas<Indices, Sim_Trait>::self_pointer PolytropicIdealGas<Indices, Sim_Trait>::_instance = NULL;

template <typename Indices, typename Sim_Trait>
typename PolytropicIdealGas<Indices, Sim_Trait>::self_reference PolytropicIdealGas<Indices, Sim_Trait>::Instance(void)
{
    if (_instance != NULL)
        return *_instance;
    else
    {
        _instance = new self_type;
        return *_instance;
    }
}

template <typename Indices, typename Sim_Trait>
PolytropicIdealGas<Indices, Sim_Trait>::PolytropicIdealGas(void)
{}

template <typename Indices, typename Sim_Trait>
PolytropicIdealGas<Indices, Sim_Trait>::~PolytropicIdealGas(void)
{}

template <typename Indices, typename Sim_Trait>
void PolytropicIdealGas<Indices, Sim_Trait>::Pressure(void)
{
#ifndef NDEBUG
    matrix_reference Data(mem_manager_type::Instance().Data);
    matrix_reference GData(mem_manager_type::Instance().GData);

    for (size_t i = 0; i < Data.size1(); ++i ) {
	    assert( Data(i, E) > 0 );
	    assert( Data(i, RHO) > 0 );
    }
    
    for (size_t i = 0; i < GData.size1(); ++i ) {
	    assert( GData(i, E) > 0 );
	    assert( GData(i, RHO) > 0 );
    }
	    
#endif
    
    gamma = mem_manager_type::Instance().LoadParameter("GAMMA");
    if (gamma != gamma) {
	    gamma = 1.4;
    }
    mem_manager_type::Instance().SaveParameter("GAMMA", gamma, true);

    const value_type const_gamma = gamma;

    p = (const_gamma - 1.) * num::element_prod(rho, e);
    gp = (const_gamma - 1.) * num::element_prod(grho, ge);
    
}

}

#endif
