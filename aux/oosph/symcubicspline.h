#ifndef OOSPHSYMCUBICSPLINE_H
#define OOSPHSYMCUBICSPLINE_H

#include <cstdlib>
#include <boost/mpl/at.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <cmath>

#include "kernel.h"
#include "manager.h"

namespace oosph
{

namespace mpl = boost::mpl;
namespace num = boost::numeric::ublas;

/**
 * \brief Cubic Spline Kernel with averaged 
 * 	  smoothing length
 * 
 * \author Pascal Bauer pbauer@phim.unibe.ch
 * 	   Andreas Reufer andreas.reufer@space.unibe.ch
 */
template <typename Indices, typename Sim_Trait>
class SymCubicSpline : public Kernel<SymCubicSpline<Indices, Sim_Trait>, Sim_Trait>
{
public:

    enum Index { X = mpl::at_c<Indices, 0>::type::value,
                 Y = X + 1,
                 Z = X + 2,
                 H = mpl::at_c<Indices, 1>::type::value};

    typedef SymCubicSpline  self_type;
    typedef SymCubicSpline& self_reference;
    typedef SymCubicSpline* self_pointer;

    typedef Kernel<SymCubicSpline<Indices, Sim_Trait>, Sim_Trait> parent_type;
    typedef Kernel<SymCubicSpline<Indices, Sim_Trait>, Sim_Trait>& parent_reference;
    typedef Kernel<SymCubicSpline<Indices, Sim_Trait>, Sim_Trait>* parent_pointer;

    typedef Sim_Trait sim_trait;

    typedef Manager<sim_trait>  manager_type;
    typedef Manager<sim_trait>& manager_reference;
    typedef Manager<sim_trait>* manager_pointer;

    typedef typename parent_type::value_type value_type;
    typedef typename parent_type::vector_type vector_type;
    typedef typename parent_type::zero_vector zero_vector;
    typedef typename parent_type::matrix_row matrix_row;
    typedef typename parent_type::matrix_row_range matrix_row_range;
    typedef typename parent_type::const_matrix_row_range const_matrix_row_range;
    
    static self_reference Instance(void);
    static value_type K(const matrix_row& i, const matrix_row& j);
    static vector_type Grad(const matrix_row& i, const matrix_row& j);

protected:

    SymCubicSpline(void);
    ~SymCubicSpline(void);

private:

    static self_pointer _instance;
    static manager_reference Manager;

    static value_type r, h, q, w;
    static vector_type rvector;

    static const num::range Pos;
};

template <typename Indices, typename Sim_Trait>
typename SymCubicSpline<Indices, Sim_Trait>::value_type
SymCubicSpline<Indices, Sim_Trait>::r = 0.;

template <typename Indices, typename Sim_Trait>
typename SymCubicSpline<Indices, Sim_Trait>::value_type
SymCubicSpline<Indices, Sim_Trait>::h = 0.;

template <typename Indices, typename Sim_Trait>
typename SymCubicSpline<Indices, Sim_Trait>::value_type
SymCubicSpline<Indices, Sim_Trait>::q = 0.;

template <typename Indices, typename Sim_Trait>
typename SymCubicSpline<Indices, Sim_Trait>::value_type
SymCubicSpline<Indices, Sim_Trait>::w = 0.;

template <typename Indices, typename Sim_Trait>
typename SymCubicSpline<Indices, Sim_Trait>::vector_type
SymCubicSpline<Indices, Sim_Trait>::rvector = zero_vector(3);

template <typename Indices, typename Sim_Trait>
const num::range SymCubicSpline<Indices, Sim_Trait>::Pos(num::range(X, X + 3) );

template <typename Indices, typename Sim_Trait>
typename SymCubicSpline<Indices, Sim_Trait>::manager_reference SymCubicSpline<Indices, Sim_Trait>::Manager(manager_type::Instance() );

template <typename Indices, typename Sim_Trait>
typename SymCubicSpline<Indices, Sim_Trait>::self_pointer SymCubicSpline<Indices, Sim_Trait>::_instance = NULL;

template <typename Indices, typename Sim_Trait>
inline typename SymCubicSpline<Indices, Sim_Trait>::self_reference SymCubicSpline<Indices, Sim_Trait>::Instance(void)
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
SymCubicSpline<Indices, Sim_Trait>::SymCubicSpline(void)
{}

template <typename Indices, typename Sim_Trait>
SymCubicSpline<Indices, Sim_Trait>::~SymCubicSpline(void)
{}

template <typename Indices, typename Sim_Trait>
inline typename SymCubicSpline<Indices, Sim_Trait>::value_type SymCubicSpline<Indices, Sim_Trait>::K(const matrix_row& i, const matrix_row& j)
{
    assert( i(H) >= 0. );
    assert( j(H) >= 0. );

    h = 0.5 * ( i(H) + j(H) );

    r = sqrt(	  (j(X) - i(X))*(j(X) - i(X))
		+ (j(Y) - i(Y))*(j(Y) - i(Y))
		+ (j(Z) - i(Z))*(j(Z) - i(Z)) );

    q = r / h;

    if (1. <= q && q <= 2.) {
	    return  ( 1. / ( M_PI*h*h*h ) ) * (0.25*(2. - q)*(2. - q)*(2. - q) );
    } else if (q < 1.) {
	    return  ( 1. / ( M_PI*h*h*h ) ) * (1. - 1.5*q*q + 0.75*q*q*q);
    } else {
	    return 0.;
    }

}

template <typename Indices, typename Sim_Trait>
inline typename SymCubicSpline<Indices, Sim_Trait>::vector_type SymCubicSpline<Indices, Sim_Trait>::Grad(const matrix_row& i, const matrix_row& j)
{
    assert(i.index() < i.data().size1() );
    assert(j.index() < j.data().size1() );

    assert( i(H) >= 0. );
    assert( j(H) >= 0. );
    
    const const_matrix_row_range pos1(i, Pos);
    const const_matrix_row_range pos2(j, Pos);
    
    h = 0.5 * ( i(H) + j(H) );

    rvector[0] = j(X) - i(X);
    rvector[1] = j(Y) - i(Y);
    rvector[2] = j(Z) - i(Z);
    
    r = sqrt(		  rvector[0]*rvector[0]
			+ rvector[1]*rvector[1]
			+ rvector[2]*rvector[2] );

    if (r <= 0.) {
	    return 0. * rvector;
    }

    q = r / h;

    if (1. <= q && q <= 2. ) {
	    w = - (3. / (4.*M_PI*h*h*h*h) ) * (2. - q) * (2. - q);
    } else if (q < 1.) {
	    w =  ( 1. / ( M_PI*h*h*h*h ) ) * ( -3.*q  + 2.25*q*q);
    } else {
	    w = 0.;
    }

    return - w * rvector / r;
}

}
;

#endif
