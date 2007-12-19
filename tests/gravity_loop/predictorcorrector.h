#ifndef OOSPHPREDICTORCORRECTOR_H
#define OOSPHPREDICTORCORRECTOR_H

#include <cstdlib>

#include <boost/mpl/at.hpp>
#include <boost/mpl/vector_c.hpp>

#include "integrator.h"

namespace oosph
{

namespace mpl = boost::mpl;

template <typename Indices, typename SimTrait>
class FOPredictorCorrector : public Integrator<FOPredictorCorrector<Indices, SimTrait>, SimTrait>
{
public:

    enum Index{X  = mpl::at_c<Indices, 0>::type::value,
               VX  = mpl::at_c<Indices, 1>::type::value,
               OX  = mpl::at_c<Indices, 2>::type::value,
               OVX = mpl::at_c<Indices, 3>::type::value,
               PX  = mpl::at_c<Indices, 4>::type::value,
               PVX = mpl::at_c<Indices, 5>::type::value };

    typedef FOPredictorCorrector  self_type;
    typedef FOPredictorCorrector& self_reference;
    typedef FOPredictorCorrector* self_pointer;

    typedef Integrator<self_type, SimTrait>  parent_type;
    typedef Integrator<self_type, SimTrait>& parent_reference;
    typedef Integrator<self_type, SimTrait>* parent_pointer;

    typedef typename parent_type::manager_type manager_type;
    typedef typename parent_type::manager_reference manager_reference;
    typedef typename parent_type::manager_pointer manager_pointer;

    typedef typename manager_type::mem_manager_type mem_manager_type;
    typedef typename manager_type::mem_manager_reference mem_manager_reference;
    typedef typename manager_type::mem_manager_pointer mem_manager_pointer;

    typedef typename SimTrait::matrix_type   matrix_type;
    typedef typename SimTrait::matrix_column matrix_column;

    typedef typename SimTrait::value_type value_type;

    FOPredictorCorrector(void);
    ~FOPredictorCorrector(void);

    void BootStrap(void);
    void Predictor(const value_type dt);
    void Corrector(const value_type dt);


protected:
private:

    static mem_manager_reference Mem;

    matrix_column x;
    matrix_column vx;
    matrix_column ox;
    matrix_column ovx;
    matrix_column px;
    matrix_column pvx;

};

template <typename Indices, typename SimTrait>
typename FOPredictorCorrector<Indices, SimTrait>::mem_manager_reference
FOPredictorCorrector<Indices, SimTrait>::Mem(mem_manager_type::Instance() );

template <typename Indices, typename SimTrait>
FOPredictorCorrector<Indices, SimTrait>::FOPredictorCorrector(void) :
        x(Mem.Data, X) ,
        vx(Mem.Data, VX),
        ox(Mem.Data, OX),
        ovx(Mem.Data, OVX),
        px(Mem.Data, PX),
        pvx(Mem.Data, PVX)
{}

template <typename Indices, typename SimTrait>
FOPredictorCorrector<Indices, SimTrait>::~FOPredictorCorrector(void)
{}

template <typename Indices, typename SimTrait>
void FOPredictorCorrector<Indices, SimTrait>::BootStrap(void)
{
    pvx = vx;
    ovx = vx;
}

template <typename Indices, typename SimTrait>
void FOPredictorCorrector<Indices, SimTrait>::Predictor(const value_type dt)
{
    px  = x + 0.5 * dt *(3. * vx - ovx);
    ox  = x;
    ovx = vx;
    x   = px;
}

template <typename Indices, typename SimTrait>
void FOPredictorCorrector<Indices, SimTrait>::Corrector(const value_type dt)
{
    x = ox + 0.5 * dt * (vx + ovx);
}

//*** Vector Predictor Corector for first order ODE
template <typename Indices, typename SimTrait>
class FOVecPredictorCorrector : Integrator<FOVecPredictorCorrector<Indices, SimTrait>, SimTrait>
{
public:
    enum Index{X  = mpl::at_c<Indices, 0>::type::value,
               Y  = X + 1,
               Z  = X + 2,
               VX  = mpl::at_c<Indices, 1>::type::value,
               VY  = VX + 1,
               VZ  = VX + 2,
               OX  = mpl::at_c<Indices, 2>::type::value,
               OY  = OX + 1,
               OZ  = OX + 2,
               OVX = mpl::at_c<Indices, 3>::type::value,
               OVY = OVX + 1,
               OVZ = OVX + 2,
               PX  = mpl::at_c<Indices, 4>::type::value,
               PY  = PX + 1,
               PZ  = PX + 2,
               PVX = mpl::at_c<Indices, 5>::type::value,
               PVY = PVX + 1,
               PVZ = PVX + 2};

    typedef typename SimTrait::value_type value_type;

    void BootStrap(void)
    {
        xint.BootStrap();
        yint.BootStrap();
        zint.BootStrap();
    }

    void Predictor(const value_type dt)
    {
        xint.Predictor(dt);
        yint.Predictor(dt);
        zint.Predictor(dt);
    }

    void Corrector(const value_type dt)
    {
        xint.Corrector(dt);
        yint.Corrector(dt);
        zint.Correcotr(dt);
    }

protected:

private:

    typedef mpl::vector_c<size_t, X, VX, OX, OVX, PX, PVX> XInd;
    typedef mpl::vector_c<size_t, Y, VY, OY, OVY, PY, PVY> YInd;
    typedef mpl::vector_c<size_t, Z, VZ, OZ, OVZ, PZ, PVZ> ZInd;

    typedef FOPredictorCorrector<XInd, SimTrait> XInt;
    typedef FOPredictorCorrector<YInd, SimTrait> YInt;
    typedef FOPredictorCorrector<ZInd, SimTrait> ZInt;

    XInt xint;
    YInt yint;
    ZInt zint;

};


template <typename Indices, typename SimTrait>
class SOPredictorCorrector : public Integrator<SOPredictorCorrector<Indices, SimTrait>, SimTrait>
{
public:

    enum Index {X  = mpl::at_c<Indices, 0>::type::value,
                VX  = mpl::at_c<Indices, 1>::type::value,
                AX  = mpl::at_c<Indices, 2>::type::value,
                OX  = mpl::at_c<Indices, 3>::type::value,
                OVX = mpl::at_c<Indices, 4>::type::value,
                OAX = mpl::at_c<Indices, 5>::type::value,
                PX  = mpl::at_c<Indices, 6>::type::value,
                PVX = mpl::at_c<Indices, 7>::type::value,
                PAX = mpl::at_c<Indices, 8>::type::value};

    typedef SOPredictorCorrector self_type;
    typedef SOPredictorCorrector& self_reference;
    typedef SOPredictorCorrector* self_pointer;

    typedef Integrator<self_type, SimTrait>  parent_type;
    typedef Integrator<self_type, SimTrait>& parent_reference;
    typedef Integrator<self_type, SimTrait>* parent_pointer;

    typedef typename parent_type::manager_type manager_type;
    typedef typename parent_type::manager_reference manager_reference;
    typedef typename parent_type::manager_pointer manager_pointer;

    typedef typename manager_type::mem_manager_type mem_manager_type;
    typedef typename manager_type::mem_manager_reference mem_manager_reference;
    typedef typename manager_type::mem_manager_pointer mem_manager_pointer;

    typedef typename SimTrait::matrix_type   matrix_type;
    typedef typename SimTrait::matrix_column matrix_column;

    typedef typename SimTrait::value_type value_type;

    SOPredictorCorrector(void);
    ~SOPredictorCorrector(void);

    void BootStrap(void);
    void Predictor(value_type dt);
    void Corrector(value_type dt);

protected:

private:

    static mem_manager_reference Mem;

    matrix_column x;
    matrix_column vx;
    matrix_column ax;
    matrix_column ox;
    matrix_column ovx;
    matrix_column oax;
    matrix_column px;
    matrix_column pvx;
    matrix_column pax;
};

template <typename Indices, typename SimTrait>
typename SOPredictorCorrector<Indices, SimTrait>::mem_manager_reference
SOPredictorCorrector<Indices, SimTrait>::Mem(mem_manager_type::Instance() );

template <typename Indices, typename SimTrait>
SOPredictorCorrector<Indices, SimTrait>::SOPredictorCorrector(void) :
        x(Mem.Data, X),
        vx(Mem.Data, VX),
        ax(Mem.Data, AX),
        ox(Mem.Data, OX),
        ovx(Mem.Data, OVX),
        oax(Mem.Data, OAX),
        px(Mem.Data, PX),
        pvx(Mem.Data, PVX),
        pax(Mem.Data, PAX)
{
    // Empty Constructor
}

template <typename Indices, typename SimTrait>
SOPredictorCorrector<Indices, SimTrait>::~SOPredictorCorrector(void)
{
    // Empty Destructor
}

template <typename Indices, typename SimTrait>
void SOPredictorCorrector<Indices, SimTrait>::BootStrap(void)
{
    pax = ax;
    pvx = vx;
    oax = ax;
    ovx = vx;
}

template <typename Indices, typename SimTrait>
void SOPredictorCorrector<Indices, SimTrait>::Predictor(value_type dt)
{
    pvx = vx + 0.5 * dt * (3. * ax - oax);
    px  = x + 0.5 * dt * (3. * vx - ovx);
    ox  = x;
    ovx = vx;
    oax = ax;
    vx  = pvx;
    x   = px;
}

template <typename Indices, typename SimTrait>
void SOPredictorCorrector<Indices, SimTrait>::Corrector(value_type dt)
{
    vx = ovx + 0.5 * dt * (ax + oax);
    x  = ox  + 0.5 * dt * (pvx +ovx);
}

//*** Vector Predictor Corrector
template <typename Indices, typename SimTrait>
class SOVecPredictorCorrector : Integrator<SOVecPredictorCorrector<Indices, SimTrait>, SimTrait>
{

    enum Index {X = mpl::at_c<Indices, 0>::type::value,
                Y = X + 1,
                Z = X + 2,
                VX = mpl::at_c<Indices, 1>::type::value,
                VY = VX + 1,
                VZ = VX + 2,
                AX = mpl::at_c<Indices, 2>::type::value,
                AY = AX + 1,
                AZ = AX + 2,
                OX = mpl::at_c<Indices, 3>::type::value,
                OY = OX + 1,
                OZ = OX + 2,
                OVX = mpl::at_c<Indices, 4>::type::value,
                OVY = OVX + 1,
                OVZ = OVX + 2,
                OAX = mpl::at_c<Indices, 5>::type::value,
                OAY = OAX + 1,
                OAZ = OAX + 2,
                PX = mpl::at_c<Indices, 6>::type::value,
                PY = PX + 1,
                PZ = PX + 2,
                PVX = mpl::at_c<Indices, 7>::type::value,
                PVY = PVX + 1,
                PVZ = PVX + 2,
                PAX = mpl::at_c<Indices, 8>::type::value,
                PAY = PAX + 1,
                PAZ = PAX + 2};

public:

    typedef typename SimTrait::value_type value_type;

    void BootStrap(void)
    {
        xint.BootStrap();
        yint.BootStrap();
        zint.BootStrap();
    }

    void Predictor(const value_type dt)
    {
        xint.Predictor(dt);
        yint.Predictor(dt);
        zint.Predictor(dt);
    }

    void Corrector(const value_type dt)
    {
        xint.Corrector(dt);
        yint.Corrector(dt);
        zint.Corrector(dt);
    }

protected:
private:

    typedef mpl::vector_c<size_t, X, VX, AX, OX, OVX, OAX, PX, PVX, PAX> XInd;
    typedef mpl::vector_c<size_t, Y, VY, AY, OY, OVY, OAY, PY, PVY, PAY> YInd;
    typedef mpl::vector_c<size_t, Z, VZ, AZ, OZ, OVZ, OAZ, PZ, PVZ, PAZ> ZInd;

    typedef SOPredictorCorrector<XInd, SimTrait> XInt;
    typedef SOPredictorCorrector<YInd, SimTrait> YInt;
    typedef SOPredictorCorrector<ZInd, SimTrait> ZInt;

    XInt xint;
    YInt yint;
    ZInt zint;
};

}

#endif
