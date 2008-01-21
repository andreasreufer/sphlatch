#ifndef OOSPH_RANKSPACE_H
#define OOSPH_RANKSPACE_H

#include <cmath>
#include <cstdlib>

#include <vector>
#include <algorithm>

#include <mpi.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/mpl/at.hpp>

#include "memorymanager.h"
#include "neighbours.h"
#include "bucket.h"
#include "costzone.h"

namespace oosph
{

namespace num = boost::numeric::ublas;

template <typename Indices, typename SimTrait>
class RankSpace : public Neighbours<RankSpace<Indices, SimTrait>, SimTrait>
{
public:

    template <typename T1, typename T2>
    friend class ParticleFunction;

    enum Index { X = mpl::at_c<Indices, 0>::type::value,
                 Y = X + 1,
                 Z = X + 2,
                 H = mpl::at_c<Indices, 1>::type::value,
                 NONEIGH = mpl::at_c<Indices, 2>::type::value };

    typedef RankSpace self_type;
    typedef RankSpace& self_reference;
    typedef RankSpace* self_pointer;

    typedef typename SimTrait::value_type value_type;

    typedef typename SimTrait::particle_container_type  particle_container_type;
    typedef typename SimTrait::particle_container_reference particle_container_reference;
    typedef typename SimTrait::particle_container_pointer particle_container_pointer;

    typedef Bucket<size_t> rankspace_bucket_type;
    typedef rankspace_bucket_type::cell_type cell_type;

    typedef typename SimTrait::matrix_reference matrix_reference;
    typedef typename SimTrait::matrix_row matrix_row;
    typedef typename SimTrait::matrix_row_range matrix_row_range;
    typedef typename SimTrait::const_matrix_row_range const_matrix_row_range;

    typedef typename SimTrait::matrix_column matrix_column;

    typedef typename SimTrait::converter_type converter_type;

    typedef num::matrix<size_t> index_matrix_type;
    typedef num::matrix_row<index_matrix_type> index_matrix_row;

    typedef MemoryManager<SimTrait> memory_type;
    
    typedef mpl::vector_c<size_t, X> CostZoneIndex;
    typedef CostZone<CostZoneIndex, SimTrait> costzone_type;

    class Ordering;

    static self_reference Instance(void);

    self_reference operator()(void);
    particle_container_reference operator[](const size_t index);
    static particle_container_reference GetNeighbours(const size_t index);

    static size_t ppc;

    costzone_type& CostZone;

protected:

    RankSpace(void);
    ~RankSpace(void);

    static particle_container_type neighbours;
private:

    static rankspace_bucket_type Bucket;

    static self_pointer _instance;
    static matrix_reference Data;
    static matrix_reference GData;

    static index_matrix_type MinMax;
    static size_t last;
    static std::vector<std::vector<size_t> > Transcription;
    static std::vector<std::vector<size_t> > Sorted;

};

template <typename Indices, typename SimTrait>
typename RankSpace<Indices, SimTrait>::self_pointer
RankSpace<Indices, SimTrait>::_instance = NULL;

template <typename Indices, typename SimTrait>
typename RankSpace<Indices, SimTrait>::matrix_reference
RankSpace<Indices, SimTrait>::Data(memory_type::Instance().Data);

template <typename Indices, typename SimTrait>
size_t RankSpace<Indices, SimTrait>::ppc = 2;

template <typename Indices, typename SimTrait>
typename RankSpace<Indices, SimTrait>::rankspace_bucket_type
RankSpace<Indices, SimTrait>::Bucket;

template <typename Indices, typename SimTrait>
typename RankSpace<Indices, SimTrait>::matrix_reference
RankSpace<Indices, SimTrait>::GData(memory_type::Instance().GData);

template <typename Indices, typename SimTrait>
typename RankSpace<Indices, SimTrait>::index_matrix_type
RankSpace<Indices, SimTrait>::MinMax;

template <typename Indices, typename SimTrait>
size_t RankSpace<Indices, SimTrait>::last = static_cast<size_t>(-1);

template <typename Indices, typename SimTrait>
typename RankSpace<Indices, SimTrait>::particle_container_type
RankSpace<Indices, SimTrait>::neighbours;

template <typename Indices, typename SimTrait>
std::vector<std::vector<size_t> >
RankSpace<Indices, SimTrait>::Transcription(3, std::vector<size_t>() );

template <typename Indices, typename SimTrait>
std::vector<std::vector<size_t> >
RankSpace<Indices, SimTrait>::Sorted(3, std::vector<size_t>() );

template <typename Indices, typename SimTrait>
typename RankSpace<Indices, SimTrait>::self_reference
RankSpace<Indices, SimTrait>::Instance(void)
{
    if (_instance != NULL)
        return *_instance;
    else
    {
        _instance = new RankSpace;
        return *_instance;
    }
}

template <typename Indices, typename SimTrait>
RankSpace<Indices, SimTrait>::RankSpace(void) : CostZone(costzone_type::Instance())
{
    ppc = 2;
}

template <typename Indices, typename SimTrait>
RankSpace<Indices, SimTrait>::~RankSpace(void)
{}


template <typename Indices, typename SimTrait>
typename RankSpace<Indices, SimTrait>::self_reference
RankSpace<Indices, SimTrait>::operator()(void)
{
    const size_t PSize = Data.size1();
    const size_t GSize = GData.size1();
    const size_t TSize = PSize + GSize;
    static const size_t SIZE = MPI::COMM_WORLD.Get_size();

    const size_t dim = 3;

    /** Start measuring time */

    CostZone.StartTimer();

    /** Initialize the Index Vectors */

    {
        std::vector<size_t>& Current(Sorted[0]);
        Current.resize(TSize);
        for (size_t i = 0; i < TSize; ++i)
        {
            Current[i] = i;
        }
        Sorted[1] = Current;
        Sorted[2] = Current;
    }

    /** Sort the Index Vectors */

    const size_t DimIdx[3] =
        {
            X, Y, Z
        };

    for (size_t i = 0; i < dim; ++i)
    {
        std::sort(Sorted[i].begin(), Sorted[i].end(), Ordering(DimIdx[i]) );
    }

    /** Transcription of the Index Vectors */

    for (size_t i = 0; i < Transcription.size(); ++i)
    {
        Transcription[i].resize(TSize);
    }

    for (size_t i = 0; i < dim; ++i)
    {
        std::vector<size_t>& CT(Transcription[i]);
        const std::vector<size_t>& CS(Sorted[i]);

        for (size_t j = 0; j < TSize; ++j)
        {
            CT[CS[j]] = j;
        }
    }

    /** Create Bucket and fill up the Values */

    const size_t side = converter_type::convert(cbrt( static_cast<value_type>(TSize) / static_cast<value_type>(ppc) ) );
    //const size_t face = side * side;
    const size_t cells = side * side * side;

    const value_type v_side = static_cast<value_type>(side);
    const value_type v_TSize = static_cast<value_type>(TSize);

    const value_type factor = v_side / ( (1. + (0.001/SIZE) ) * v_TSize );

    Bucket.clear();
    Bucket.resize(cells);

    const size_t idx_help[] =
        {
            1, side, side * side
        };

    for (size_t i = 0; i < TSize; ++i)
    {
        size_t index = 0;

        for (size_t j = 0; j < dim; ++j)
        {
            assert(factor * static_cast<value_type>(Transcription[j][i]) < side);
            index += converter_type::convert(factor * static_cast<value_type>(Transcription[j][i]) ) * idx_help[j];
        }

        assert (index < cells);

        Bucket[index].push_back(i);
    }

    /** Prearrange for Neighbour Search */

    MinMax.resize(PSize, 6);

    /** Do the Minimums first */

    const matrix_column PCol[] =
        {
            matrix_column(Data, X),
            matrix_column(Data, Y),
            matrix_column(Data, Z)
        };

    const matrix_column GCol[] =
        {
            matrix_column(GData, X),
            matrix_column(GData, Y),
            matrix_column(GData, Z)
        };


    const size_t MaxIndex = TSize - 1;

    size_t LMM[6];

    for (size_t i = 0; i < PSize; ++i)
    {
        const value_type diff = 2 * Data(i, H);
        const value_type MM[] =
            {
                Data(i, X) - diff,
                Data(i, Y) - diff,
                Data(i, Z) - diff,
                Data(i, X) + diff,
                Data(i, Y) + diff,
                Data(i, Z) + diff
            };

        for (size_t dim = 0; dim < 3; ++dim)
        {

            const std::vector<size_t>& ST(Sorted[dim]);
            const std::vector<size_t>& CT(Transcription[dim]);

            // Start Performance Improvement

            {
                size_t start = 0;
                size_t end = CT[i];

                while ( (end - start) > 1 )
                {
                    const size_t median = start + (end - start) / 2;
                    const size_t& I = ST[median];
                    const bool ghost = (I >= PSize);
                    const size_t d_idx = ghost ? I - PSize : I;
                    const matrix_column& CC = ghost ? GCol[dim] : PCol[dim];
                    const value_type& value = CC(d_idx);
                    const bool lower = (value < MM[dim]);
                    if (lower)
                    {
                        start = median;
                    }
                    else
                    {
                        end = median;
                    }

                    //                     std::cout << "start: " << start << "\tEnd: " << end << "\tdiff: " << (end - start) << std::endl;
                }

                //                 std::cout << "Start / End: " << ST[start] << "\t" << ST[end] << std::endl;
                //                 std::cout << start << "\t" << end << std::endl;
                LMM[dim] = start;

            }

            // End Performance Improvement

            //             {
            //                 bool cont = true;
            //                 size_t& min = LMM[dim] = CT[i];
            //
            //                 while (cont)
            //                 {
            //                     const size_t& I = ST[min];
            //                     const bool ghost = (I >= PSize);
            //                     const size_t d_idx = ghost ? I - PSize : I;
            //                     const matrix_column& CC = ghost ? GCol[dim] : PCol[dim];
            //                     const double& value = CC(d_idx);
            //                     const bool mina = (value < MM[dim] );
            //                     const bool minb = (min == 0);
            //
            //                     if (!(mina || minb) )
            //                         --min;
            //                     else
            //                     {
            //                         if (mina)
            //                             ++min;
            //                         cont = false;
            //                     }
            //                 }
            //             }

            // Start Performance Improvement
            {
                size_t start = CT[i];
                size_t end   = MaxIndex;

                while ( (end - start) > 1 )
                {
                    const size_t median = start + (end - start) / 2;
                    const size_t& I = ST[median];
                    const bool ghost = (I >= PSize);
                    const size_t d_idx = ghost ? I - PSize : I;
                    const matrix_column& CC = ghost ? GCol[dim] : PCol[dim];
                    const value_type& value = CC(d_idx);
                    const bool higher = (value > MM[dim + 3]);
                    if (higher)
                    {
                        end = median;
                    }
                    else
                    {
                        start = median;
                    }
                }
                
                LMM[dim + 3] = start;

            }
            // End Performance Improvement

        }

        const size_t dim_idx[] =
            {
                0, 1, 2, 0, 1, 2
            };

        for (size_t dim = 0; dim < 6; ++dim)
        {
            MinMax(i, dim) = Sorted[dim_idx[dim]][LMM[dim]];
//             std::cout << "MinMax:\t" << MinMax(i, dim) << std::endl;
        }

    }

    /** Check the Buckets ... */

    //     for (size_t i = 0; i < cells; ++i)
    //     {
    //         cell_type& Current = Bucket[i];
    //
    //         for (cell_type::const_iterator P = Current.begin(); P != Current.end(); ++P)
    //         {
    //             const size_t& idx = *P;
    //             const size_t cidx = idx < PSize ? idx : idx - PSize;
    //             const matrix_reference CData = idx < PSize ? Data : GData;
    //             matrix_row CRow(CData, cidx);
    //             std::cout << matrix_row_range(CRow, num::range(X, X + 3) ) << std::endl;
    //         }
    //         std::cout << "-------------------------------" << std::endl;
    //     }

    /** Check the Bucket */

    //     for (size_t i = 0; i < cells; ++i)
    //     {
    //         cell_type& Current = Bucket[i];
    //
    //         for (cell_type::const_iterator P = Current.begin(); P != Current.end(); ++P)
    //         {
    //
    //             const size_t& idx = *P;
    //             const size_t cidx = idx < PSize ? idx : idx - PSize;
    //             const matrix_reference CData = idx < PSize ? Data : GData;
    //             matrix_row PRow(CData, cidx);
    //             matrix_row_range PPos(PRow, num::range(X, X + 3) );
    //             std::cout << "H: " << PRow(H) << std::endl;
    //
    //             for (cell_type::const_iterator Q = Current.begin(); Q != Current.end(); ++Q)
    //             {
    //                 const size_t& idx = *Q;
    //                 const size_t cidx = idx < PSize ? idx : idx - PSize;
    //                 const matrix_reference CData = idx < PSize ? Data : GData;
    //                 matrix_row QRow(CData, cidx);
    //                 matrix_row_range QPos(QRow, num::range(X, X + 3) );
    // 		const double d = num::norm_2(QPos - PPos);
    //                 std::cout << (d > PRow(H) ? "*\t" : "\t" ) << num::norm_2(QPos - PPos) << std::endl;
    //             }
    //         }
    //         std::cout << "-----------------------------" << std::endl;
    //     }

    /** Cechk Transcription and Sort */

    //     for (size_t i = 0; i < TSize; ++i)
    //     {
    //         for (size_t j = 0; j < 3; ++j)
    //         {
    //             std::cout << i << "\t" << Transcription[j][i] << "\t" << Sorted[j][i] << "\t" << Transcription[j][Sorted[j][i]] << "\t" << Sorted[j][Transcription[j][i]] << std::endl;
    //         }
    //     }

    /** Stop measuring time */

    CostZone.StopTimer();
    
    return *this;
}

template <typename Indices, typename SimTrait>
typename RankSpace<Indices, SimTrait>::particle_container_reference
RankSpace<Indices, SimTrait>::GetNeighbours(const size_t index)
{
    if (last == index)
        return neighbours;
    else
        last = index;

    neighbours.clear();

    static matrix_reference Data (memory_type::Instance().Data);
    static matrix_reference GData (memory_type::Instance().GData);

    assert (index < Data.size1() );

    const matrix_row Current(Data, index);
    const const_matrix_row_range Pos(Current, num::range(X, X + 3) );
    const value_type dist = 2 * Current(H);

    const size_t PSize = Data.size1();
    const size_t GSize = GData.size1();
    const size_t TSize = PSize + GSize;
    static const size_t SIZE = MPI::COMM_WORLD.Get_size();

    const size_t side = converter_type::convert(cbrt( static_cast<value_type>(TSize) / static_cast<value_type>(ppc) ) );
    const size_t face = side * side;
    ////////////////const size_t cells = side * side * side;

    const value_type factor = static_cast<value_type>(side) / ( ( 1. + (0.001/SIZE) ) * static_cast<value_type>(TSize) );

    //     #define CHECKNEIGHBOURS
#ifndef CHECKNEIGHBOURS

    const size_t bminmax[] =
        {
	      (size_t) ( /*(value_type) */Transcription[0][MinMax(index, 0)] ) * factor,
	      (size_t) ( /*(value_type)*/ Transcription[1][MinMax(index, 1)] ) * factor,
	      (size_t) ( /*(value_type)*/ Transcription[2][MinMax(index, 2)] ) * factor,
	      (size_t) ( /*(value_type)*/ Transcription[0][MinMax(index, 3)] ) * factor,
	      (size_t) ( /*(value_type)*/ Transcription[1][MinMax(index, 4)] ) * factor,
	      (size_t) ( /*(value_type)*/ Transcription[2][MinMax(index, 5)] ) * factor
        };
#endif
#ifdef CHECKNEIGHBOURS

    const size_t bminmax[] =
        {
            0, 0, 0, side - 1, side -1, side - 1
        };
#endif

    for (size_t k = bminmax[0]; k <= bminmax[3]; ++k)
    {
        for (size_t j = bminmax[1]; j <= bminmax[4]; ++j)
        {
            for (size_t i = bminmax[2]; i <= bminmax[5]; ++i)
            {
                const size_t idx = k + j * side + i * face;
                const cell_type& Current(Bucket[idx]);

                for (cell_type::const_iterator pos = Current.begin(); pos != Current.end(); ++pos)
                {
                    const size_t cindex = *pos < PSize ? *pos : *pos - PSize;
                    const matrix_reference CData = *pos < PSize ? Data : GData;
                    const matrix_row CRow(CData, cindex);
                    const const_matrix_row_range CPos(CRow, num::range(X, X + 3) );

                    const value_type dd = sqrt( (Pos(0) - CPos(0) ) * (Pos(0) - CPos(0) ) + (Pos(1) - CPos(1) ) * (Pos(1) - CPos(1) ) + (Pos(2) - CPos(2) ) * (Pos(2) - CPos(2) ) );

                    if (dd <= dist)
                    {
                        //                         if (Data(index, 0) == CData(cindex, 0) )
                        //                             continue;

                        neighbours.push_back(matrix_row(CData, cindex) );
                    }
                }
            }
        }
    }
    
    // Account number of neighbours for better smoothing length estimation
    // Add one, because current particle also counts as neighbour
    Data(index, NONEIGH) = neighbours.size() + 1;
    return neighbours;
}

template <typename Indices, typename SimTrait>
inline typename RankSpace<Indices, SimTrait>::particle_container_reference
RankSpace<Indices, SimTrait>::operator[](const size_t index)
{
    return GetNeighbours(index);
}

template <typename Indices, typename SimTrait>
class RankSpace<Indices, SimTrait>::Ordering
{
public:

    typedef typename SimTrait::matrix_column matrix_column;
    typedef typename SimTrait::matrix_row    matrix_row;

    Ordering(const size_t index) : N(index), PSize(Data.size1()) , GSize(GData.size1() )
    {}
    ;

    const bool operator()(const size_t& a, const size_t& b)
    {
        const matrix_reference Data_1 = a >= PSize ? GData : Data;
        const matrix_reference Data_2 = b >= PSize ? GData : Data;

        const size_t index_1 = a >= PSize ? a - PSize : a;
        const size_t index_2 = b >= PSize ? b - PSize : b;

        return Data_1(index_1, N) < Data_2(index_2, N);
    }

protected:

private:
    const size_t PSize;
    const size_t GSize;
    const size_t N;
};

}
;

#endif
