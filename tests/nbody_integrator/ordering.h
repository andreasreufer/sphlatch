#ifndef OOSPH_ORDERING_H
#define OOSPH_ORDERING_H

#include <cmath>
#include <iostream>
#include <bitset>

#include <boost/dynamic_bitset.hpp>

namespace oosph
{

struct OrderingTag
    {}
;

struct Cartesian    : OrderingTag
    {}
;
struct Morton       : OrderingTag
    {}
;
struct PeanoHilbert : OrderingTag
    {}
;

template <int N>
struct OrderingInfo
{
    static const size_t depth;
    static const size_t side;
    static const size_t face;
    static const size_t cells;
};

template <int N>
const size_t OrderingInfo<N>::depth = N;

template <int N>
const size_t OrderingInfo<N>::side = lrint(pow(2, N) );

template <int N>
const size_t OrderingInfo<N>::face = side * side;

template <int N>
const size_t OrderingInfo<N>::cells = side * side * side;

template <typename TypeTag, int N>
struct OrderingDispatch
{
    OrderingDispatch(void)
    {}

    static size_t Index(const size_t& i, const size_t& j, const size_t& k)
    {
        return 0;
    }

    static void Coordinate(const size_t& Index, size_t& i, size_t& j, size_t& k)
    {
        i = 0;
        j = 0;
        k = 0;
    }

    static OrderingInfo<N> Cube;
}
;

template <int N>
struct OrderingDispatch<Cartesian, N>
{
    typedef OrderingInfo<N> cube;

    OrderingDispatch(void)
    {}

    static size_t Index(const size_t& i, const size_t& j, const size_t& k)
    {
        return i + j * cube::side + k * cube::side * cube::side;
    }

    static void Coordinate(const size_t& Index, size_t& i, size_t& j, size_t& k)
    {

        k = Index - (Index % (cube::face) );
        j = Index - k - (Index % cube::side);
        i = Index - k - j;

        j /= cube::side;
        k /= cube::face;

        return;
    }
}
;

template <int N>
struct OrderingDispatch<Morton, N>
{
    OrderingDispatch(void)
    {}

    static size_t Index(const size_t& i, const size_t& j, const size_t& k)
    {
        boost::dynamic_bitset<> IndexBits(3 * N);

        std::vector<boost::dynamic_bitset<> > Bit;

        Bit.push_back(boost::dynamic_bitset<>(N, i) );
        Bit.push_back(boost::dynamic_bitset<>(N, j) );
        Bit.push_back(boost::dynamic_bitset<>(N, k) );

        for (size_t b = 0; b < N; ++b)
        {
            for (size_t a = 0; a < 3; ++a)
            {
                IndexBits[3 * b + a] = Bit[a][b];
            }
        }

        return static_cast<size_t>(IndexBits.to_ulong() );
    }

    static void Coordinate(const size_t& Index, size_t& i, size_t& j, size_t& k)
    {

        boost::dynamic_bitset<> IndexBits(3 * N, Index);

        std::vector<boost::dynamic_bitset<> > Bit(3, boost::dynamic_bitset<>(N) );

        for (size_t b = 0; b < N; ++b)
        {
            for (size_t a= 0; a < 3; ++a)
            {
                Bit[a][b] = IndexBits[3 * b + a];
            }
        }

        i = static_cast<int>(Bit[0].to_ulong() );
        j = static_cast<int>(Bit[1].to_ulong() );
        k = static_cast<int>(Bit[2].to_ulong() );
    }

    static OrderingInfo<N> Cube;
};

template <int N>
struct OrderingDispatch<PeanoHilbert, N>
{
    OrderingDispatch(void)
    {}

    static size_t Index(const size_t& i, const size_t& j, const size_t& k)
    {

        return 3;
    }
    static void Coordinate(const size_t& Index, size_t& i, size_t& j, size_t& k)
    {
        boost::dynamic_bitset<> IndexBits(3 * N, Index);
        std::bitset<3> Q;

        int x, y, z;
        int X, Y, Z;

        X = Y = Z = 0;
        x = y = z = 0;

        for (size_t a = 0; a < N; ++a)
        {
            for (size_t b = 0; b < 3; ++b)
            {
                Q[b] = IndexBits[a * 3 + b];
            }

            x = X;
            y = Y;
            z = Z;

            const size_t type = Q.to_ulong();

            switch (type)
            {
            case 0:
                X = z;
                Y = y;
                Z = x;
                break;
            case 1:
                X = y;
                Y = x;
                Z = z + static_cast<int>( pow( 2, a ) );
                break;
            case 2:
                X = y;
                Y = x + static_cast<int>( pow( 2, a ) );
                Z = z + static_cast<int>( pow( 2, a ) );
                break;
            case 3:
                X = x;
                Y = static_cast<int>( pow( 2, a + 1 ) ) - 1 - y;
                Z = static_cast<int>( pow( 2, a ) ) - 1 - z;
                break;
            case 4:
                X = x + static_cast<int>( pow( 2, a ) );
                Y = static_cast<int>( pow( 2, a + 1 ) ) - 1 - y;
                Z = static_cast<int>( pow( 2, a ) ) - 1 - z;
                break;
            case 5:
                X = static_cast<int>( pow( 2, a + 1 ) ) - 1 - y;
                Y = static_cast<int>( pow( 2, a + 1 ) ) - 1 - x;
                Z = z + static_cast<int>( pow( 2, a ) );
                break;
            case 6:
                X = static_cast<int>( pow( 2, a + 1 ) ) - 1 - y;
                Y = static_cast<int>( pow( 2, a ) ) - 1 - x;
                Z = z + static_cast<int>( pow( 2, a ) );
                break;
            case 7:
                X = static_cast<int>( pow( 2, a + 1 ) ) - 1 - z;
                Y = y;
                Z = static_cast<int>( pow( 2, a ) ) - 1 - x;
                break;

            }
        }

        i = static_cast<size_t>(X);
        j = static_cast<size_t>(Y);
        k = static_cast<size_t>(Z);

        return;
    }

    static OrderingInfo<N> Cube;
};

// ***** Ordering Class ******
template <typename TypeTag, int N>
class Ordering
{
public:

    const size_t operator()(const size_t i, const size_t j, const size_t k);
    const void   operator()(const size_t Index, size_t& i, size_t& j, size_t& k);

protected:
private:

    OrderingDispatch<TypeTag, N> Dispatch;

    static const size_t depth;
    static const size_t side;
    static const size_t cells;

};

template <typename TypeTag, int N>
const size_t Ordering<TypeTag, N>::depth = N;

template <typename TypeTag, int N>
const size_t Ordering<TypeTag, N>::side = lrint(pow(2, N));

template <typename TypeTag, int N>
const size_t Ordering<TypeTag, N>::cells = side * side * side;

template <typename TypeTag, int N>
inline const size_t Ordering<TypeTag, N>::operator()(const size_t i, const size_t j, const size_t k)
{
    return OrderingDispatch<TypeTag, N>::Index(i, j, k);
}

template <typename TypeTag, int N>
inline const void Ordering<TypeTag, N>::operator()(const size_t Index, size_t& i, size_t& j, size_t& k)
{
    return OrderingDispatch<TypeTag, N>::Coordinate(Index, i, j, k);
}


}
;

#endif
