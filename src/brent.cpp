#ifndef SPHLATCH_BRENT
#define SPHLATCH_BRENT

#include "typedefs.h"

namespace sphlatch {
template<typename _funcT>
class BrentRooter {
public:
   BrentRooter() :
      maxItr(50),
      eps(3.e-8),
      nan(std::numeric_limits<fType>::quiet_NaN()),
      quiet(false)
   { }
   ~BrentRooter() { }

private:
   class noConvergence {
public:
      noConvergence(const fType _a, const fType _fa,
                    const fType _b, const fType _fb,
                    const fType _c, const fType _fc)
      {
         std::cerr << "no convergence in BrentRooter:\n"
                   << "   a: " << _a << "   fa: " << _fa << "\n"
                   << "   b: " << _b << "   fb: " << _fb << "\n"
                   << "   c: " << _c << "   fc: " << _fc << "\n";
      }
   };
   class rootNotBracketed {
public:
      rootNotBracketed(const fType _a, const fType _fa,
                       const fType _b, const fType _fb)
      {
         std::cerr << "not bracketed  in BrentRooter:\n"
                   << "   a: " << _a << "   fa: " << _fa << "\n"
                   << "   b: " << _b << "   fb: " << _fb << "\n";
      }
   };

public:
   fType operator()(const fType _x1, const fType _x2, const fType _tol)
   {
      fType a, b, c = 0., d = 0., e = 0.;

      a = _x1;
      b = _x2;

      fType fa, fb, fc, P, Q, R, S, min1, min2, tol1, xm;
      fa = f(a);
      fb = f(b);
      if (((fa > 0.) && (fb > 0.)) or ((fa < 0.) && (fb < 0.)))
      {
         if (quiet)
            return(nan);
         else
            throw rootNotBracketed(a, fa, b, fb);
      }

      fc = fb;

      for (size_t itr = 1; itr <= maxItr; itr++)
      {
         if (((fb > 0.) && (fc > 0.)) or ((fb < 0.) && (fc < 0.)))
         {
            c  = a;
            fc = fa;
            e  = d = b - a;
         }
         if (fabs(fc) < fabs(fb))
         {
            a  = b;
            b  = c;
            c  = a;
            fa = fb;
            fb = fc;
            fc = fa;
         }

         tol1 = 2.* eps* fabs(b) + 0.5 * _tol;
         xm   = 0.5 * (c - b);

         if (fabs(xm) <= tol1 or fb == 0.0)
         {
            f(b);
            return(b);
         }

         if ((fabs(e) >= tol1) && (fabs(fa) > fabs(fb)))
         {
            S = fb / fa;
            if (a == c)
            {
               P = 2.0 * xm * S;
               Q = 1. - S;
            }
            else
            {
               Q = fa / fc;
               R = fb / fc;
               P = S * (2. * xm * Q * (Q - R) - (b - a) * (R - 1.));
               Q = (Q - 1.) * (R - 1.) * (S - 1.);
            }
            if (P > 0.)
               Q = -Q;
            P    = fabs(P);
            min1 = 3. * xm * Q - fabs(tol1 * Q);
            min2 = fabs(e * Q);
            if (2. * P < (min1 < min2 ? min1 : min2))
            {
               e = d;
               d = P / Q;
            }
            else
            {
               d = xm;
               e = d;
            }
         }
         else
         {
            d = xm;
            e = d;
         }
         a  = b;
         fa = fb;
         if (fabs(d) > tol1)
            b += d;
         else
            b += (xm > 0. ? fabs(tol1) : -fabs(tol1));
         fb = f(b);
      }
      if (not quiet)
         throw noConvergence(a, fa, b, fb, c, fc);
      else
         return(nan);
   }

   _funcT f;

private:
   size_t maxItr;
   fType  eps, nan;
public:
   bool quiet;
};
}

#endif
