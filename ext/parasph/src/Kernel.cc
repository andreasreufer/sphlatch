/****************3*************************************************************
 * ParaSPH -- Version 24.11.2003                                              *
 *----------------------------------------------------------------------------*
 * File:      Kernel.cc                                                       *
 * Purpose:   This class computes the beta spline kernel W(v) for values      *
 *            v = 0 .. 2. You can access its values directly using the f-     *
 *            methods, or the get-methods to get an interpolated value (which *
 *            runs faster on some machines). The latter case requires v*v as  *
 *            input. If you need some details look at appendix B of my thesis.*
 *****************************************************************************/
#ifndef KERNEL_CC
#define KERNEL_CC

#include <assert.h>
#include <math.h>
#include <iomanip>

//#include "Def.cc"

typedef double ktype;
//typedef ftype ktype;

class Kernel {
private:
  static const ktype A, B, C, D, E, F, G, H, I, J, K, L, M, N;
  enum { wd = 40001 };
  ktype dw;

  ktype W[wd+1], gW[wd+1], fm[wd+1];

public:
  static const ktype vmax, v2max;

  // fm and gW will always be multiplicated by r_vec/r. In order to save
  // calculation of r, we do it here. We just have to multiply the returned
  // value by 1/h.
  ktype ffm(const ktype &v) const {
    assert(v <= vmax);
    ktype v2 = v*v, v3 = v*v*v, v4 = v*v3, v5 = v*v4;
    if (v >= 1.0) return -J / v + K * v2 - L * v3 + M * v4 - N * v5;
    else return H * v2 - G * v4 + I * v5;
  }

  ktype fW(const ktype &v) const {
    assert(v <= vmax);
    if (v >= 1.0) return D*(2.-v)*(2.-v)*(2.-v);
    return A + v*v*(C*v - B);
  }
  
  ktype fgW(const ktype &v) const {
    assert(v <= vmax);
    if (v >= 1.0) return F - F/v - C*v;
    return E*v - F;
  }
  
  Kernel() {
    ktype v;
    dw = wd * 0.25;
    for (int i = 0; i <= wd; i++) {
      v = sqrt(i / dw);
      fm[i] = ffm(v);
      W[i]  = fW(v);
      gW[i] = fgW(v);
    }
  }

  ktype getFm(const ktype &v2) const {
    int index = (int)(v2*dw);
    // new assertion: array +1 element (since v2 maybe 4 + epsilon)
    assert(index <= wd);
    return fm[index] + (fm[index+1] - fm[index]) * (v2 * dw - index);
  }

  ktype getW(const ktype &v2) const {
    int index = (int)(v2*dw);
    // new assertion: array +1 element (since v2 maybe 4 + epsilon)
    assert(index <= wd);
    return W[index] + (W[index+1] - W[index]) * (v2 * dw - index);
  }

  ktype getgW(const ktype &v2) const {
    int index = (int)(v2*dw);
    // new assertion: array +1 element (since v2 maybe 4 + epsilon)
    assert(index <= wd);    
    return gW[index] + (gW[index+1] - gW[index]) * (v2 * dw - index); 
  }
};

const ktype Kernel::A  = 1./ 3.14159265358979323846;
const ktype Kernel::B  = A * 1.5;
const ktype Kernel::C  = A * 0.75;
const ktype Kernel::D  = A * 0.25;
const ktype Kernel::E  = A * 2.25;
const ktype Kernel::F  = A * 3.0;

const ktype Kernel::G  = 4. / 3.;
const ktype Kernel::H  = 1.2;
const ktype Kernel::I  = 0.5;
const ktype Kernel::J  = 1. / 15.;
const ktype Kernel::K  = 8. / 3.;
const ktype Kernel::L  = 3.;
const ktype Kernel::M  = 1.2;
const ktype Kernel::N  = 1. / 6.;

const ktype Kernel::vmax  = 2.0;
const ktype Kernel::v2max = vmax*vmax;

#endif
