#ifndef OOSPH_COSTZONE_H
#define OOSPH_COSTZONE_H

#include <iostream>
#include <cmath>
#include <limits>
#include <vector>
#include <deque>

#include <boost/mpl/at.hpp>

#include "loadbalance.h"
#include "domainmap.h"
#include "bucket.h"
#include "memorymanager.h"
#include "communicationmanager.h"

#include "ordering.h"

namespace oosph
{

namespace mpl = boost::mpl;

template <typename Indices, typename SimTrait>
class CostZone : public LoadBalance<CostZone<Indices, SimTrait>, SimTrait>
{
public:

    enum Index {X = mpl::at_c<Indices, 0>::type::value,
                Y = X + 1,
                Z = X + 2};

    enum config {default_depth = 4};

    typedef CostZone self_type;
    typedef CostZone& self_reference;
    typedef CostZone* self_pointer;

    typedef LoadBalance<CostZone<Indices, SimTrait>, SimTrait> parent_type;
    typedef parent_type& parent_reference;
    typedef parent_type* parent_pointer;

    typedef typename SimTrait::matrix_type matrix_type;
    typedef typename SimTrait::matrix_reference matrix_reference;
    typedef typename SimTrait::matrix_pointer matrix_pointer;
    
    typedef typename SimTrait::matrix_row matrix_row;
    typedef typename SimTrait::matrix_column matrix_column;

    typedef typename SimTrait::value_type value_type;
    typedef typename SimTrait::value_reference value_reference;
    typedef typename SimTrait::value_pointer value_pointer;

    typedef typename SimTrait::vector_type vector_type;
    typedef typename SimTrait::vector_reference vector_reference;
    typedef typename SimTrait::vector_pointer vector_pointer;
    
    typedef typename SimTrait::index_vector_type index_vector_type;
    typedef typename SimTrait::index_vector_reference index_vector_reference;
    typedef typename SimTrait::index_vector_pointer index_vector_pointer;

    typedef typename SimTrait::index_list_type index_list_type;
    typedef typename SimTrait::index_list_reference index_list_reference;
    typedef typename SimTrait::index_list_pointer index_list_pointer;

    typedef MemoryManager<SimTrait> memory_manager_type;
    typedef MemoryManager<SimTrait>& memory_manager_reference;
    typedef MemoryManager<SimTrait>* memory_manager_pointer;

    typedef CommunicationManager<SimTrait> communication_manager_type;
    typedef CommunicationManager<SimTrait>& communication_manager_reference;
    typedef CommunicationManager<SimTrait>* communication_manager_pointer;

    typedef std::vector<index_vector_type> domain_index_vector_type;
    typedef std::vector<index_vector_type>& domain_index_vector_reference;
    typedef std::vector<index_vector_type>* domain_index_vector_pointer;

    typedef std::vector<index_list_type> domain_index_list_type;
    typedef std::vector<index_list_type>& domain_index_list_reference;
    typedef std::vector<index_list_type>* domain_index_list_pointer;

    typedef Bucket<size_t> bucket_type;
    typedef Bucket<size_t>& bucket_reference;
    typedef Bucket<size_t>* bucket_pointer;

    typedef std::numeric_limits<value_type> limits;

    typedef Ordering<Cartesian, default_depth> cartesian_ordering_type;
    typedef Ordering<PeanoHilbert, default_depth> peano_hilbert_ordering_type;

    typedef DomainMap domain_map_type;
    typedef DomainMap& domain_map_reference;
    typedef DomainMap* domain_map_pointer;

    typedef typename SimTrait::converter_type converter_type;


protected:


private:

    static self_pointer _instance;
    static size_t depth;
    static size_t side;
    static size_t cells;

    static domain_map_type domain_map;
    static bucket_type particle_bucket;

    static memory_manager_reference memory_manager;
    static communication_manager_reference communication_manager;
    static cartesian_ordering_type cartesian_ordering;
    static peano_hilbert_ordering_type peano_hilbert_ordering;



public:

    static self_reference Instance(void);

    domain_index_vector_reference CreateDomainIndexVector(void);
    domain_index_vector_reference CreateDomainGhostIndexVector(void);

    vector_type GetLowerCorner();
    vector_type GetUpperCorner();
    value_type  GetCellSize();

    void StartTimer(void);
    void StopTimer(void);

    value_type LocalPart(void);

    static domain_index_vector_reference DomainGhostIndexVector;

    void resize(const size_t size);


protected:

    CostZone(const size_t _depth = default_depth);
    ~CostZone(void);


private:

    bool CreateParticleBucket(bucket_reference Bucket);
    bool CreateDomainMap(domain_map_reference map);

    double TimerStartTime;
    double TimerTotalTime;

    size_t TotalParticles;

    bool TimerRunning;

    std::deque<vector_type> LoadBalanceFactor;

};

template <typename Indices, typename SimTrait>
typename CostZone<Indices, SimTrait>::domain_index_vector_reference
CostZone<Indices, SimTrait>::DomainGhostIndexVector(parent_type::gdiv);

template <typename Indices, typename SimTrait>
typename CostZone<Indices, SimTrait>::cartesian_ordering_type
CostZone<Indices, SimTrait>::cartesian_ordering;

template <typename Indices, typename SimTrait>
typename CostZone<Indices, SimTrait>::peano_hilbert_ordering_type
CostZone<Indices, SimTrait>::peano_hilbert_ordering;

template <typename Indices, typename SimTrait>
typename CostZone<Indices, SimTrait>::domain_map_type
CostZone<Indices, SimTrait>::domain_map;

template <typename Indices, typename SimTrait>
typename CostZone<Indices, SimTrait>::bucket_type
CostZone<Indices, SimTrait>::particle_bucket;

template <typename Indices, typename SimTrait>
typename CostZone<Indices, SimTrait>::memory_manager_reference
CostZone<Indices, SimTrait>::memory_manager(memory_manager_type::Instance() );

template <typename Indices, typename SimTrait>
typename CostZone<Indices, SimTrait>::communication_manager_reference
CostZone<Indices, SimTrait>::communication_manager(communication_manager_type::Instance() );

template <typename Indices, typename SimTrait>
typename CostZone<Indices, SimTrait>::self_pointer
CostZone<Indices, SimTrait>::_instance = NULL;

template <typename Indices, typename SimTrait>
typename CostZone<Indices, SimTrait>::self_reference
CostZone<Indices, SimTrait>::Instance(void)
{
    if (_instance != NULL)
        return *_instance;
    else
    {
        _instance = new CostZone;
        return *_instance;
    }
}

template <typename Indices, typename SimTrait>
size_t CostZone<Indices, SimTrait>::depth = default_depth;

template <typename Indices, typename SimTrait>
size_t CostZone<Indices, SimTrait>::side = lrint(pow(2, default_depth));

template <typename Indices, typename SimTrait>
size_t CostZone<Indices, SimTrait>::cells = side * side * side;

template <typename Indices, typename SimTrait>
CostZone<Indices, SimTrait>::CostZone(const size_t _depth)
{
    const size_t SIZE = MPI::COMM_WORLD.Get_size();
    
    depth = _depth;
    side = lrint(pow(2, depth));
    cells = side * side * side;
    domain_map.resize(cells);
    particle_bucket.resize(cells);

    TimerRunning = false;
    TimerStartTime = 0;

    // Initalize with != 0 or initializing will run havoc
    TimerTotalTime = 1;

    // Determine how far back Loadbalance should average out
    const size_t average_interval = 4;

    vector_type BalancedVector;
    BalancedVector.resize(SIZE);
    for (size_t i = 0; i < SIZE; i++) {
	    BalancedVector[i] = 1./SIZE;
    }

    LoadBalanceFactor.resize(average_interval, BalancedVector);
}

template <typename Indices, typename SimTrait>
CostZone<Indices, SimTrait>::~CostZone(void)
{}

template <typename Indices, typename SimTrait>
typename CostZone<Indices, SimTrait>::domain_index_vector_reference
CostZone<Indices, SimTrait>::CreateDomainIndexVector(void)
{
    const size_t SIZE = MPI::COMM_WORLD.Get_size();
    
    CreateParticleBucket();
    CreateDomainMap();

    parent_type::div.resize( SIZE );

    for (size_t i = 0; i < SIZE; ++i)
    {
        parent_type::div[i].clear();
    }

    for (size_t i = 0; i < cells; ++i)
    {
        index_list_reference S(particle_bucket[i]);
        index_vector_reference D(parent_type::div[domain_map[i]]);
        D.insert(D.end(), S.begin(), S.end() );
        S.clear();
    }
        
    return parent_type::div;
};

template <typename Indices, typename SimTrait>
typename CostZone<Indices, SimTrait>::domain_index_vector_reference
CostZone<Indices, SimTrait>::CreateDomainGhostIndexVector(void)
{
    const size_t SIZE = MPI::COMM_WORLD.Get_size();
    const size_t RANK = MPI::COMM_WORLD.Get_rank();

    CreateParticleBucket();

    const size_t face = side * side;

    domain_index_vector_reference DIVR(parent_type::gdiv);

    DIVR.resize( SIZE );

    for (size_t i = 0; i < SIZE; ++i)
    {
        DIVR[i].clear();
    }

    assert( domain_map.size() == cells);

    for (size_t z = 0; z < side; ++z)
    {
        for (size_t y = 0; y < side; ++y )
        {
            for (size_t x = 0; x < side; ++x)
            {
                const size_t index = cartesian_ordering(x, y, z);

                assert( index < cells);

                if (domain_map[index] != RANK )
                    continue;

                for (int k = -1; k <= 1; ++k)
                {
                    for (int j = -1; j <= 1; ++j)
                    {
                        for (int i = -1; i <= 1; ++i)
                        {
                            const int a = x + i;
                            const int b = y + j;
                            const int c = z + k;

                            if (a < 0 || b < 0 || c < 0)
                                continue;
                            if (a >= (int)side || b >= (int)side || c >= (int)side)
                                continue;
                            if (a == (int)x && b == (int)y && c == (int)z)
                                continue;

                            const size_t nindex = a + b * side + c * face;
                            const size_t cindex = x + y * side + z * face;

                            assert(nindex <= cells);
                            assert(cindex <= cells);

                            if (domain_map[nindex] == RANK )
                                continue;

                            index_list_reference S(particle_bucket[cindex]);
                            index_vector_reference D(DIVR[domain_map[nindex]]);

                            D.insert(D.end(), S.begin(), S.end() );
                        }
                    }
                }
            }
        }
    }

    for (size_t i = 0; i < SIZE; ++i)
    {
        std::sort(DIVR[i].begin(), DIVR[i].end() );
        DIVR[i].resize(std::unique(DIVR[i].begin(), DIVR[i].end() ) - DIVR[i].begin() );
    }

    return parent_type::gdiv;
}

template <typename Indices, typename SimTrait>
void CostZone<Indices, SimTrait>::resize(const size_t size)
{
    depth = size;
    side = lrint(pow(2, depth) );
    cells = side * side * side;
    domain_map.clear();
    domain_map.resize(cells);
    particle_bucket.resize(cells);
}

template <typename Indices, typename SimTrait>
bool CostZone<Indices, SimTrait>::CreateParticleBucket(bucket_reference Bucket = particle_bucket)
{
    matrix_reference Data(memory_manager.Data);

    Bucket.clear();
    Bucket.resize(cells);

    for (size_t i = 0; i < cells; ++i)
    {
        Bucket[i].clear();
    }

    const size_t DSIZE = Data.size1();

    /* Change min/max */
    const value_type min_x = communication_manager.template Min<X>();
    const value_type max_x = communication_manager.template Max<X>();
    const value_type min_y = communication_manager.template Min<Y>();
    const value_type max_y = communication_manager.template Max<Y>();
    const value_type min_z = communication_manager.template Min<Z>();
    const value_type max_z = communication_manager.template Max<Z>();

    // This is broken for unsymmetric problems:
    // the interval may be made up from a difference of different dimensions!
    const value_type min = min_x < ( ( min_y < min_z ) ? min_y : min_z ) ? min_x : ( ( min_y < min_z ) ? min_y : min_z ) ;

    const value_type max = max_x > ( ( max_y > max_z ) ? max_y : max_z ) ? max_x : ( ( max_y > max_z ) ? max_y : max_z ) ;

    const value_type Interval = fabs(max - min);
    const value_type NInt = static_cast<value_type>(side) / Interval;
    const value_type NIntMin = min * NInt;

    const int max_index = (int)side - 1;

    for (size_t cp = 0; cp < DSIZE; ++cp)
    {

        assert(DSIZE <= Data.size1() );

        matrix_row Current(Data, cp);

        const size_t i = std::max( 0, std::min( max_index, converter_type::convert(NInt * Current(X) - NIntMin) ) );
        const size_t j = std::max( 0, std::min( max_index, converter_type::convert(NInt * Current(Y) - NIntMin) ) );
        const size_t k = std::max( 0, std::min( max_index, converter_type::convert(NInt * Current(Z) - NIntMin) ) );

        particle_bucket[cartesian_ordering(i, j, k)].push_back(cp);
    }

    return true;
}


template <typename Indices, typename SimTrait>
bool CostZone<Indices, SimTrait>::CreateDomainMap(domain_map_reference map = domain_map)
{
    const size_t SIZE = MPI::COMM_WORLD.Get_size();
    
    const size_t BuffSize = 4096; // Must be smaller than least number of cells!
    
    TotalParticles = 0;
   
    /* Let's determine the cell sizes */
    std::vector<size_t> GlobalCellsizes(cells);

    /* We'll have to this the buffered way, as too big arrays tend to segfault */
    for (size_t offset = 0; offset < cells; offset += BuffSize) {
	    
	    // Fill the buffer
	    size_t *LocalCellsizesBuff = new size_t[BuffSize];
	    size_t *GlobalCellsizesBuff = new size_t[BuffSize];

	    for (size_t i = 0; i < BuffSize; i++) {
		    LocalCellsizesBuff[i] = particle_bucket[i + offset].size();
		    GlobalCellsizesBuff[i] = 0;
	    }

	    /* Reduce */
	    MPI::COMM_WORLD.Reduce(LocalCellsizesBuff, GlobalCellsizesBuff, BuffSize, MPI::UNSIGNED_LONG, MPI::SUM, 0);
	    MPI::COMM_WORLD.Bcast(GlobalCellsizesBuff, BuffSize, MPI::UNSIGNED_LONG, 0);

	    /* Write back and count particles */
	    for (size_t i = 0; i < BuffSize; i++) {
		    GlobalCellsizes[i + offset] = GlobalCellsizesBuff[i];
		    TotalParticles += GlobalCellsizesBuff[i];
	    }

	    delete[] LocalCellsizesBuff;
	    delete[] GlobalCellsizesBuff;
    }

    /* Exchange usage information about last time step */

    double GlobalTotalTimes[SIZE];
    MPI::COMM_WORLD.Allgather(&TimerTotalTime, 1, MPI::DOUBLE, &GlobalTotalTimes, 1, MPI::DOUBLE);

    /* Reset timer */
    TimerTotalTime = 0;

    vector_type LoadBalanceLatestFactor;
    LoadBalanceLatestFactor.resize(SIZE);
    for (size_t i = 0; i < SIZE; i++) {
	    LoadBalanceLatestFactor[i] = 1. / GlobalTotalTimes[i];
	    /** OOSPH_NOLOADBALANCE disables loadbalancing */

#ifdef OOSPH_NOLOADBALANCE 
	    LoadBalanceLatestFactor[i] = 1.;
#endif

    }

    LoadBalanceLatestFactor = LoadBalanceLatestFactor / norm_1(LoadBalanceLatestFactor);

    LoadBalanceFactor.pop_front();
    LoadBalanceFactor.push_back(LoadBalanceLatestFactor);
    
    /* Now determine domain sizes */

    vector_type LoadBalanceCurrentFactor;
    std::vector<size_t>  DomainSize;
    
    LoadBalanceCurrentFactor.resize(SIZE);
    DomainSize.resize(SIZE);

    for (size_t i = 0; i < SIZE; i++) {
	    LoadBalanceCurrentFactor[i] = 0;
	    for (size_t j = 0; j < LoadBalanceFactor.size(); j++) {
		    LoadBalanceCurrentFactor[i] += ( LoadBalanceFactor[j] )[i];
	    }
	    LoadBalanceCurrentFactor[i] = LoadBalanceCurrentFactor[i] / LoadBalanceFactor.size();
	    DomainSize[i] = lrint( LoadBalanceCurrentFactor[i] * TotalParticles );
    }
    
    /* Calculate domain map */

    size_t CurDomain = 0;
    size_t MaxDomain = SIZE - 1;

    size_t CurParticles = 0;
    size_t DistParticles = 0;

    domain_map.resize(cells);

    for (size_t i = 0; i < cells; ++i) {
	    size_t a, b, c;
	    peano_hilbert_ordering(i, a, b, c);
	    const size_t cart_i = cartesian_ordering(a, b, c);

	    // If ( addition of current cell completes distribution OR
	    //  domain is saturated )
	    //  AND current domain is not the last one 
	    //  AND current domain has more than 0 particles
	    // then begin filling next domain

	    if ( ( DistParticles + GlobalCellsizes[cart_i] == TotalParticles || CurParticles >= DomainSize[CurDomain] )
			   && CurDomain < MaxDomain && CurParticles > 0 ) {
		    CurParticles = 0;
		    ++CurDomain;
	    }

	    CurParticles += GlobalCellsizes[cart_i];
	    DistParticles += GlobalCellsizes[cart_i];
	    
	    domain_map[cart_i] = CurDomain;
    }

    return true;
}


template <typename Indices, typename SimTrait>
void CostZone<Indices, SimTrait>::StartTimer(void)
{
	if (!TimerRunning) {
		TimerStartTime = MPI_Wtime();
		TimerRunning = true;
	}
}

template <typename Indices, typename SimTrait>
void CostZone<Indices, SimTrait>::StopTimer(void)
{
	value_type TimerStopTime;
		
	if (TimerRunning) {
		TimerStopTime	= MPI_Wtime();
		TimerRunning = false;
		
		TimerTotalTime += ( TimerStopTime - TimerStartTime );
	}
}

template <typename Indices, typename SimTrait>
typename CostZone<Indices, SimTrait>::value_type
CostZone<Indices, SimTrait>::LocalPart(void)
{
	matrix_reference Data(memory_manager.Data);
	return ( (value_type)Data.size1() / TotalParticles );
}

/*template <typename Indices, typename SimTrait>
typename CostZone<Indices, SimTrait>::value_type
CostZone<Indices, SimTrait>::LocalPart(void)
{
	matrix_reference Data(memory_manager.Data);
	return ( (value_type)Data.size1() / TotalParticles );
}*/

}
;

#endif
