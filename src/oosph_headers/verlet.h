#ifndef OOSPHVERLET_H
#define OOSPHVERLET_H

#include <cstdlib>

#include <boost/mpl/at.hpp>
#include <boost/mpl/vector_c.hpp>

#include "integrator.h"


namespace oosph
{
namespace mpl = boost::mpl;

template <typename Indices, typename SimTrait>
class Verlet : public Integrator<Verlet<Indices, SimTrait>, SimTrait>
{
public:

enum Index { X = mpl::at_c<Indices, 0>::type::value,
             VX = mpl::at_c<Indices, 1>::type::value,
             AX = mpl::at_c<Indices, 2>::type::value,
             OAX = mpl::at_c<Indices, 3>::type::value };

typedef Verlet self_type;
typedef Verlet& self_reference;
typedef Verlet* self_pointer;

typedef Integrator<self_type, SimTrait>  parent_type;
typedef Integrator<self_type, SimTrait>& parent_reference;
typedef Integrator<self_type, SimTrait>* parent_pointer;

typedef sphlatch::MemoryManager  mem_manager_type;
typedef sphlatch::MemoryManager&  mem_manager_reference;
typedef sphlatch::MemoryManager*  mem_manager_pointer;

typedef typename SimTrait::matrix_type matrix_type;
typedef typename SimTrait::matrix_column matrix_column;

typedef typename SimTrait::value_type value_type;

Verlet(void);
~Verlet(void);

void getPos(value_type dt);
void getVel(value_type dt);

protected:

private:

static mem_manager_reference Mem;

matrix_column x;
matrix_column vx;
matrix_column ax;
matrix_column oax;
};

template <typename Indices, typename SimTrait>
typename Verlet<Indices, SimTrait>::mem_manager_reference
Verlet<Indices, SimTrait>::Mem(mem_manager_type::instance());

template <typename Indices, typename SimTrait>
Verlet<Indices, SimTrait>::Verlet(void) :
  x(Mem.Data, X),
  vx(Mem.Data, VX),
  ax(Mem.Data, AX),
  oax(Mem.Data, OAX)
{
  // Empty Constructor
}

template <typename Indices, typename SimTrait>
Verlet<Indices, SimTrait>::~Verlet(void)
{
  // Empty Destructor
}

template <typename Indices, typename SimTrait>
void Verlet<Indices, SimTrait>::getPos(value_type dt)
{
  x = x + vx*dt + ax*0.5*dt*dt;
  oax = ax;
}

template <typename Indices, typename SimTrait>
void Verlet<Indices, SimTrait>::getVel(value_type dt)
{
  vx = vx + ax*0.5*dt + oax*0.5*dt;
}

//*** Vector Predictor Corrector
template <typename Indices, typename SimTrait>
class VecVerlet : Integrator<VecVerlet<Indices, SimTrait>, SimTrait>
{
enum Index { X = mpl::at_c<Indices, 0>::type::value,
             Y = X + 1,
             Z = X + 2,
             VX = mpl::at_c<Indices, 1>::type::value,
             VY = VX + 1,
             VZ = VX + 2,
             AX = mpl::at_c<Indices, 2>::type::value,
             AY = AX + 1,
             AZ = AX + 2,
             OAX = mpl::at_c<Indices, 3>::type::value,
             OAY = OAX + 1,
             OAZ = OAX + 2 };

public:

typedef typename SimTrait::value_type value_type;

void getPos(const value_type dt)
{
  xint.getPos(dt);
  yint.getPos(dt);
  zint.getPos(dt);
}

void getVel(const value_type dt)
{
  xint.getVel(dt);
  yint.getVel(dt);
  zint.getVel(dt);
}

private:

typedef mpl::vector_c<size_t, X, VX, AX, OAX> XInd;
typedef mpl::vector_c<size_t, Y, VY, AY, OAY> YInd;
typedef mpl::vector_c<size_t, Z, VZ, AZ, OAZ> ZInd;

typedef Verlet<XInd, SimTrait> XInt;
typedef Verlet<YInd, SimTrait> YInt;
typedef Verlet<ZInd, SimTrait> ZInt;

XInt xint;
YInt yint;
ZInt zint;
};
}

#endif
