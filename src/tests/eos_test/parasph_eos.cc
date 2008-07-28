#include <math.h>

typedef double ftype;

typedef struct {
  ftype rho0, A, B, a, b, alpha, beta, u0, Eiv, Ecv;
} eosMat;

void parasphTillotson(const ftype &rho, const ftype &rhom1, const ftype &u,
                      const int &mat, ftype &P, ftype &cs, ftype &T,
                      const eosMat& m) 
{
  //eosMat m = eosTab.get(mat);
  ftype PC = 0.;
  ftype csC = 0.;
  ftype rho0m1 = 1. / m.rho0;
  ftype eta = rho * rho0m1;
  ftype mu = eta - 1.;

  ftype csmin = 0.25 * m.A * rho0m1;
  ftype Pmin = -1.e15;

  ftype c1 = u / (m.u0 * eta * eta);
  ftype c2 = 1. / (c1 + 1.);

  ///
  /// vaporization, P needed
  ///
  if (u > m.Eiv && eta < 1.)
    {
      ftype d1 = m.rho0 * rhom1;
      ftype d2 = d1 - 1.;
      ftype ex1 = exp(-m.beta * d2);
      ftype ex2 = exp(-m.alpha * d2 * d2);

      P = m.a * rho * u;
      P += ex2 * (m.b * rho * u * c2 + m.A * mu * ex1);

      cs = m.b * u * (3. * c1 + 1) * c2 * c2 + 2. * m.alpha * d2 * m.b * d1 * u * c2;
      cs += m.A * ex1 * ((2. * m.alpha * d2 + m.beta) * mu * d1 * rhom1 + rho0m1);
      cs = cs * ex2 + m.a * u;
      cs += P * rhom1 * (m.a + m.b * c2 * c2 * ex2);
      if (cs < 0.) cs = 0.;
    }

 
  ///
  /// not completely vaporized, PC needed
  /// 
  if (u < m.Ecv || eta >= 1.)
    {
      PC = (m.a + m.b * c2) * rho * u + m.A * mu + m.B * mu * mu;

      csC = m.a * u + rho0m1 * (m.A + 2. * m.B * mu) + m.b * u * (3. * c1 + 1.) * c2 * c2;
      csC += PC * rhom1 * (m.a + m.b * c2 * c2);

      ///
      /// last two conditions are redundant
      ///
      if (u > m.Eiv && u < m.Ecv && eta < 1)
        {
          ftype e1 = m.Ecv - u;
          ftype e2 = u - m.Eiv;
          ftype e3 = 1. / (m.Ecv - m.Eiv);
          P = (e2 * P + e1 * PC) * e3;
          cs = (e2 * cs + e1 * csC) * e3;
        }
      else
        {
          P = PC;
          cs = csC;
        }
    }
  if (cs < csmin) cs = csmin;
  if (P < Pmin) P = Pmin;
  cs = sqrt(cs);
}


