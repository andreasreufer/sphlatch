#ifndef OOSPH_MPI_HELP_H
#define OOSPH_MPI_HELP_H

#include <mpi.h>

namespace oosph
{

template <typename T>
struct MPI_SingleCom
{
    void min(T& t1, T& t2)
    {
        MPI::COMM_WORLD.Reduce(&t1, &t2, sizeof(T), MPI::BYTE, MPI::MIN, 0);
        MPI::COMM_WORLD.Bcast(&t1, sizeof(T), MPI::BYTE, 0);
    }
    void max(T& t1, T& t2)
    {
        MPI::COMM_WORLD.Reduce(&t1, &t2, sizeof(T), MPI::BYTE, MPI::MAX, 0);
        MPI::COMM_WORLD.Bcast(&t2, sizeof(T), MPI::BYTE, 0);
    }
}
;

template <>
struct MPI_SingleCom<float>
{
    void min(float& t1, float& t2)
    {
        MPI::COMM_WORLD.Reduce(&t1, &t2, 1, MPI::FLOAT, MPI::MIN, 0);
        MPI::COMM_WORLD.Bcast(&t2, 1, MPI::FLOAT, 0);
    }

    void max(float& t1, float& t2)
    {
        MPI::COMM_WORLD.Reduce(&t1, &t2, 1, MPI::FLOAT, MPI::MAX, 0);
        MPI::COMM_WORLD.Bcast(&t2, 1, MPI::FLOAT, 0);
    }
};

template <>
struct MPI_SingleCom<double>
{
    void min(double& t1, double& t2)
    {
        MPI::COMM_WORLD.Reduce(&t1, &t2, 1, MPI::DOUBLE, MPI::MIN, 0);
        MPI::COMM_WORLD.Bcast(&t2, 1, MPI::DOUBLE, 0);
    }

    void max(double& t1, double& t2)
    {
        MPI::COMM_WORLD.Reduce(&t1, &t2, 1, MPI::DOUBLE, MPI::MAX, 0);
        MPI::COMM_WORLD.Bcast(&t2, 1, MPI::DOUBLE, 0);
    }
};

}
;

#endif
