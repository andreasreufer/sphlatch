#ifndef OOSPHMONAGHAN_H
#define OOSPHMONAGHAN_H

#include <mpi.h>

#include <cstdlib>
#include <cmath>

#include <boost/mpl/at.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/lexical_cast.hpp>

#include "artificialviscosity.h"
#include "manager.h"


namespace oosph
{

namespace mpl = boost::mpl;
namespace num = boost::numeric::ublas;

/**
 * \brief Artificial Viscosity after Monaghan
 * 
 * \author Pascal Bauer pbauer@phim.unibe.ch
 * \author Andreas Reufer areufer@space.unibe.ch
 */

template <typename Indices, typename Sim_Trait>
class Monaghan : public ArtificialViscosity<Monaghan<Indices, Sim_Trait>, Sim_Trait>
{
public:

    enum Index { X     = mpl::at_c<Indices, 0>::type::value,
                 Y     = X + 1,
                 Z     = X + 2,
                 VX    = mpl::at_c<Indices, 1>::type::value,
                 VY    = VX + 1,
                 VZ    = VX + 2,
                 RHO   = mpl::at_c<Indices, 2>::type::value,
                 H     = mpl::at_c<Indices, 3>::type::value,
                 P     = mpl::at_c<Indices, 4>::type::value,
                 MAXAVMU = mpl::at_c<Indices, 5>::type::value};

    typedef ArtificialViscosity<Monaghan<Indices, Sim_Trait>, Sim_Trait>  parent_type;
    typedef ArtificialViscosity<Monaghan<Indices, Sim_Trait>, Sim_Trait>& parent_reference;
    typedef ArtificialViscosity<Monaghan<Indices, Sim_Trait>, Sim_Trait>* parent_pointer;

    typedef Monaghan  self_type;
    typedef Monaghan& self_reference;
    typedef Monaghan* self_pointer;

    typedef typename Sim_Trait::matrix_type memory_container_type;
    typedef typename Sim_Trait::matrix_reference memory_container_reference;
    typedef typename Sim_Trait::matrix_pointer memory_container_pointer;

    typedef Manager<Sim_Trait>  manager_type;
    typedef Manager<Sim_Trait>& manager_reference;
    typedef Manager<Sim_Trait>* manager_pointer;

    typedef typename manager_type::mem_manager_type mem_manager_type;
    typedef typename manager_type::mem_manager_reference mem_manager_reference;

    typedef typename Sim_Trait::matrix_type matrix_type;
    typedef typename Sim_Trait::matrix_reference matrix_reference;

    typedef typename parent_type::matrix_row matrix_row;
    typedef typename parent_type::value_type value_type;

    typedef typename Sim_Trait::matrix_row_range matrix_row_range;
    typedef typename Sim_Trait::matrix_column column_type;

    static self_reference Instance();

    static void PreLocal(void);
    static value_type Local(matrix_row Part1, matrix_row particle_2);
    
protected:

    Monaghan(void);
    ~Monaghan(void);

private:

    static self_pointer _instance;
    static matrix_reference Data;

    static value_type alpha;
    static value_type beta;
    static value_type gamma;
    static value_type nu;
    
    static manager_reference Manager;

};

template <typename Indices, typename Sim_Trait>
typename Monaghan<Indices, Sim_Trait>::matrix_reference
Monaghan<Indices, Sim_Trait>::Data(mem_manager_type::Instance().Data);

template <typename Indices, typename Sim_Trait>
typename Monaghan<Indices, Sim_Trait>::manager_reference Monaghan<Indices, Sim_Trait>::Manager(manager_type::Instance() );

template <typename Indices, typename Sim_Trait>
typename Monaghan<Indices, Sim_Trait>::self_pointer Monaghan<Indices, Sim_Trait>::_instance = NULL;

template <typename Indices, typename SimTrait>
typename Monaghan<Indices, SimTrait>::value_type
Monaghan<Indices, SimTrait>::alpha;

template <typename Indices, typename SimTrait>
typename Monaghan<Indices, SimTrait>::value_type
Monaghan<Indices, SimTrait>::beta;

template <typename Indices, typename SimTrait>
typename Monaghan<Indices, SimTrait>::value_type
Monaghan<Indices, SimTrait>::gamma;

template <typename Indices, typename SimTrait>
typename Monaghan<Indices, SimTrait>::value_type
Monaghan<Indices, SimTrait>::nu;

template <typename Indices, typename Sim_Trait>
typename Monaghan<Indices, Sim_Trait>::self_reference Monaghan<Indices, Sim_Trait>::Instance(void)
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
Monaghan<Indices, Sim_Trait>::Monaghan(void)
{
	mem_manager_type& MemoryManager(mem_manager_type::Instance());
	
	alpha = MemoryManager.LoadParameter("MONAGHALPHA");
	if (alpha != alpha) {
		// Setting the default value for the AV alpha
		alpha = 1.0;
	}
	MemoryManager.SaveParameter("MONAGHALPHA", alpha, true);

	beta = MemoryManager.LoadParameter("MONAGHBETA");
	if (beta != beta) {
		// Setting the default value for the AV beta
		beta = 2.0;
	}
	MemoryManager.SaveParameter("MONAGHBETA", beta, true);

	nu = MemoryManager.LoadParameter("MONAGHNU");
	if (nu != nu) {
		// Setting the default value for the AV beta
		nu = 1E-1;
	}
	MemoryManager.SaveParameter("MONAGHNU", nu, true);
	
	gamma = MemoryManager.LoadParameter("GAMMA");
	if (gamma != gamma) {
		// Setting the default value for the EOS gamma
		gamma = 1.4;
	}
	MemoryManager.SaveParameter("GAMMA", gamma, true);

}

template <typename Indices, typename Sim_Trait>
Monaghan<Indices, Sim_Trait>::~Monaghan(void)
{}

template <typename Indices, typename Sim_Trait>
void Monaghan<Indices, Sim_Trait>::PreLocal(void)
{}

template <typename Indices, typename Sim_Trait>
typename Monaghan<Indices, Sim_Trait>::value_type Monaghan<Indices, Sim_Trait>::Local(matrix_row Part1, matrix_row Part2)
{
    const value_type veldistprod =	( Part1(VX) - Part2(VX) ) * ( Part1(X) - Part2(X) )
	    			      + ( Part1(VY) - Part2(VY) ) * ( Part1(Y) - Part2(Y) )
	    			      + ( Part1(VZ) - Part2(VZ) ) * ( Part1(Z) - Part2(Z) );
		    
    assert( finite( veldistprod ) != 0 );
    
    if (veldistprod < 0.) {
        
	assert(Part1(RHO) > 0.);
        assert(Part2(RHO) > 0.);
        assert(Part1(P) > 0.);
        assert(Part2(P) > 0.);
        assert(Part1(H) > 0.);
        assert(Part2(H) > 0.);

	const value_type distsquare =   ( Part1(X) - Part2(X) ) * ( Part1(X) - Part2(X) )
				      + ( Part1(Y) - Part2(Y) ) * ( Part1(Y) - Part2(Y) )
				      + ( Part1(Z) - Part2(Z) ) * ( Part1(Z) - Part2(Z) );
	assert( distsquare >= 0. );

	const value_type h12 = ( Part1(H) + Part2(H) ) / 2.;
	const value_type rho12 = ( Part1(RHO) + Part2(RHO) ) / 2.;
	const value_type c12 = ( sqrt( gamma * Part1(P) / Part1(RHO) ) 
			 + sqrt( gamma * Part2(P) / Part2(RHO) ) ) / 2.;

	const value_type mu12 = h12 * veldistprod / ( distsquare + nu*nu*h12*h12 );
	const value_type abs_mu12 = fabs(mu12);

	assert( Part1(MAXAVMU) >= 0. );
	assert( Part2(MAXAVMU) >= 0. );
	Part1(MAXAVMU) = abs_mu12 > Part1(MAXAVMU) ? abs_mu12 : Part1(MAXAVMU);
	Part2(MAXAVMU) = abs_mu12 > Part2(MAXAVMU) ? abs_mu12 : Part2(MAXAVMU);

	return ( - alpha*mu12*c12 + beta*mu12*mu12 ) / rho12;

    } else {

        return 0.;

    }
}



/**
 * \brief Class to get the Courant Friedrichs Laax Conditions for time Integration...
 * 
 * \author Pascal Bauer pbauer@phim.unibe.ch
 * \author Andreas Reufer areufer@space.unibe.ch
 */

template <typename Indices, typename SPHTrait>
class CFL
{
public:
    enum Index {AX   = mpl::at_c<Indices, 0>::type::value,
		AY = AX + 1,
		AZ = AX + 2,
		H   = mpl::at_c<Indices, 1>::type::value,
		RHO = mpl::at_c<Indices, 2>::type::value,
		E = mpl::at_c<Indices, 3>::type::value,
		P = mpl::at_c<Indices, 4>::type::value,
		POW = mpl::at_c<Indices, 5>::type::value,
		MAXAVMU   = mpl::at_c<Indices, 6>::type::value};

    typedef CFL  self_type;
    typedef CFL& self_reference;
    typedef CFL* self_pointer;

    typedef typename SPHTrait::manager_type manager_type;
    typedef typename SPHTrait::manager_reference manager_reference;
    typedef typename SPHTrait::manager_pointer manager_pointer;

    typedef typename manager_type::mem_manager_type mem_manager_type;
    typedef typename manager_type::mem_manager_reference mem_manager_reference;
    typedef typename manager_type::mem_manager_pointer mem_manager_pointer;

    typedef typename SPHTrait::simulation_trait::matrix_row particle_type;
    typedef typename SPHTrait::simulation_trait::value_type value_type;

    static self_reference Instance(void);

    const value_type operator()(void);

protected:
    CFL(void);
    ~CFL(void);

private:
    static self_pointer _instance;
    static mem_manager_reference MemManager;
};


template <typename Indices, typename SPHTrait>
typename CFL<Indices, SPHTrait>::mem_manager_reference
CFL<Indices, SPHTrait>::MemManager(mem_manager_type::Instance() );

/** \brief Pointer to the Singleton Instance of this Class */
template <typename Indices, typename SPHTrait>
typename CFL<Indices, SPHTrait>::self_pointer CFL<Indices, SPHTrait>::_instance = NULL;

/** \brief Returns a reference to the Singleton Instance of theis Class */
template <typename Indices, typename SPHTrait>
typename CFL<Indices, SPHTrait>::self_reference CFL<Indices, SPHTrait>::Instance(void)
{
    if (_instance != NULL)
        return *_instance;
    else
    {
        _instance = new self_type;
        return *_instance;
    }
}

/** \brief Standart Constructor */
template <typename Indices, typename SPHTrait>
CFL<Indices, SPHTrait>::CFL(void)
{}

/** \brief Standart Destructor */
template <typename Indices, typename SPHTrait>
CFL<Indices, SPHTrait>::~CFL(void)
{}

/** \brief Returns the CFL-Value for the Integration Step */
template <typename Indices, typename SPHTrait>
typename SPHTrait::simulation_trait::value_type
const CFL<Indices, SPHTrait>::operator()(void)
{
    const size_t NoParticles = MemManager.Data.size1();

    value_type dt_min_acc, dt_min_av, dt_min_temp, dt_min, dt_globalmin, acc, csound;
    value_type CourantNumber, Alpha, Beta, Gamma, SpecMinEnergy;

    CourantNumber = MemManager.LoadParameter("COURANTNUMBER");
    if (CourantNumber != CourantNumber) {
	    // Setting the default value for the Courant number
	    CourantNumber = 0.5;
    }
    MemManager.SaveParameter("COURANTNUMBER", CourantNumber, true);
    
    Alpha = MemManager.LoadParameter("MONAGHALPHA");
    if (Alpha != Alpha) {
	    // Setting the default value for the AV alpha
	    Alpha = 1.;
    }
    MemManager.SaveParameter("MONAGHALPHA", Alpha, true);
    
    Beta = MemManager.LoadParameter("MONAGHBETA");
    if (Beta != Beta) {
	    // Setting the default value for the AV beta
	    Beta = 2.;
    }
    MemManager.SaveParameter("MONAGHBETA", Beta, true);
    
    Gamma = MemManager.LoadParameter("GAMMA");
    if (Gamma != Gamma) {
	    // Setting the default value for the adiabatic index gamma
	    Gamma = 1.4;
    }
    MemManager.SaveParameter("GAMMA", Gamma, true);

    SpecMinEnergy = MemManager.LoadParameter("SPECMINENERGY");
    if (SpecMinEnergy != SpecMinEnergy) {
	    // Setting the default value for the minimal
	    // specific energy ( ~ 1.3K for gamma = 5/3 and
	    // 			 ~ 0.8K for gamma = 1.4 )
	    // in the WT4 system
	    SpecMinEnergy = 1e-3;
    }
    MemManager.SaveParameter("SPECMINENERGY", SpecMinEnergy, true);

    assert (MemManager.Data.size1() > 0);
    assert (P < MemManager.Data.size2() );
    assert (RHO < MemManager.Data.size2() );

    particle_type particle0(MemManager.Data, 0);

    assert( particle0(RHO) > 0. );
    assert( particle0(P) > 0. );
    assert( particle0(MAXAVMU) >= 0. );
    assert( finite( particle0(AX) ) != 0 );
    assert( finite( particle0(AY) ) != 0 );
    assert( finite( particle0(AZ) ) != 0 );

    const value_type cAlpha = Alpha, cBeta = Beta, cGamma = Gamma, cSpecMinEnergy = SpecMinEnergy;
    acc = sqrt( particle0(AX)*particle0(AX) + particle0(AY)*particle0(AY) + particle0(AZ)*particle0(AZ) );
    assert( acc > 0. );

    csound = sqrt(cGamma * particle0(P) / particle0(RHO) );

    dt_min_acc = sqrt( particle0(H) / acc );
    assert( dt_min_acc > 0. );

    dt_min_av = particle0(H) / ( csound + 0.6*( cAlpha*csound + cBeta*particle0(MAXAVMU) ) );
    assert( dt_min_av > 0. );

    dt_min_temp = ( particle0(E) + cSpecMinEnergy ) / ( fabs( particle0(POW) ) );
    assert( dt_min_temp > 0. );

    dt_min = dt_min_acc < dt_min_av ? dt_min_acc : dt_min_av;
    dt_min = dt_min < dt_min_temp ? dt_min : dt_min_temp;

    for (size_t i = 1; i < NoParticles; ++i)
    {
        particle_type particle(MemManager.Data, i);
        
	assert( particle(RHO) > 0. );
	assert( particle(P)   > 0. );
	assert( particle(MAXAVMU) >= 0.);
	assert( finite( particle(AX) ) != 0 );
	assert( finite( particle(AY) ) != 0 );
	assert( finite( particle(AZ) ) != 0 );

	acc = sqrt( particle(AX)*particle(AX) + particle(AY)*particle(AY) + particle(AZ)*particle(AZ) );
	assert( acc > 0. );
    	csound = sqrt(cGamma * particle(P) / particle(RHO) );
    	dt_min_acc = sqrt( particle(H) / acc );
    	dt_min_av = particle(H) / ( csound + 0.6*( cAlpha*csound + cBeta*particle(MAXAVMU) ) );
	dt_min_temp = ( particle(E) + cSpecMinEnergy ) / ( fabs( particle(POW) ) );

	assert( dt_min_acc > 0. );
	assert( dt_min_av  > 0. );
	assert( dt_min_temp > 0. );

	value_type dt_min_act = dt_min_acc < dt_min_av ? dt_min_acc : dt_min_av;
	dt_min_act = dt_min_act < dt_min_temp ? dt_min_act : dt_min_temp;
	
        dt_min = dt_min_act < dt_min ? dt_min_act : dt_min;
	assert( dt_min > 0. );
    }

    MPI::COMM_WORLD.Reduce(&dt_min, &dt_globalmin, 1, MPI::DOUBLE, MPI::MIN, 0);
    MPI::COMM_WORLD.Bcast(&dt_globalmin, 1, MPI::DOUBLE, 0);

    assert( CourantNumber * dt_globalmin > 0. );
    return ( CourantNumber * dt_globalmin );

}

}

#endif
