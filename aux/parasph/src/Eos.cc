/******************************************************************************
 * ParaSPH -- Version 13.01.2004                                              *
 *----------------------------------------------------------------------------*
 * File:      Eos.cc                                                          *
 * Purpose:   Implementation of different equations of state (EOS).           *
 *            If you would like to add a new equation of state (EOS) you need *
 *            to first the new function in the format of 'eosfunc' and second *
 *            make an entry in the switch statement in init and eventually    *
 *            the 'findrho' method. Good luck :-)                             *
 * Warning:   Everything concerning ANEOS is embraced by compiler switches,   *
 *            because normal people (who do not like to use ANEOS) don't want *
 *            to link the ANEOS library to their code ...                     *
 * Changes:   If you want to add a variable (for input and/or output) in the  *
 *            call for the EOS, you need to change it at several places:      *
 *            a) The line where void (Eos::*eosfunc) is defined               *
 *            b) The 'eos' method                                             *
 *            c) The headers of all EOS subroutines                           *
 *            d) The calling subroutine (eg. in the 'findrho' method in this  *
 *               class, and in class 'Particle'                               *
 * Modifications:                                                             *
 *        -   20.12.2003 WB                                                   *
 *            Implemented aneos for arbitrary materials. Need to compile with *
 *            appropriate version of libaneos.f                               *
 *        -   13.01.2004 BN                                                   *
 *            The ANEOS input file can now be specified in the '.config' file *
 *            as well (with the variable 'input.eosDataFile'). Do not forget  *
 *            to change the variable 'simulation.eos' appropriately as well.  *
 *****************************************************************************/
#ifndef EOS_CC
#define EOS_CC

#include <math.h>

#include "Material.cc"

#ifdef ANEOS
extern "C" void aneosinit_(const char *, int);
extern "C" void aneos_(double *, double *, double *, double *, double *, 
		       double *, double *, double *, double *, double *, 
		       int *,    int *, double *, double *);
#endif

class Eos {
private:
  int              number;
  Material<eosMat> eosTab;
  void (Eos::*eosfunc)(const ftype&, const ftype &, const ftype &, const int &,
		       ftype &, ftype &, ftype &) const;

public:
  class ErrorNoNumber {
  public:
    ErrorNoNumber() {}
  };  

  // Choose eos number: 0 = ideal gas; 1 = tillotson; 2 = aneos
  void init(const int &_num, const bool &noisy, 
	    const std::string &eosDataFile) {
    number = _num;
    switch (number) {
      // add a further EOS here ... case 3: eosfunc = ........
#ifdef ANEOS
    case 2: {
      aneosinit_(eosDataFile.c_str(), eosDataFile.length());
      eosfunc = &Eos::aneos;
      if (noisy) std::cout << "Eos: Aneos" << std::endl;
      break;
    }
#endif
    case 1:
      eosfunc = &Eos::tillotson;
      eosTab.read(noisy, eosDataFile);
      if (noisy) std::cout << "Eos: Tillotson" << std::endl;
      break;
    case 0:
      eosfunc = &Eos::idealGas;
      if (noisy) std::cout << "Eos: Ideal Gas" << std::endl;
      break;
    default:
      throw ErrorNoNumber();
    }
  }

  void eos(const ftype &rho, const ftype &rhom1, const ftype &u, 
	   const int &mat, ftype &p, ftype &cs, ftype &T) const {
    (*this.*eosfunc)(rho, rhom1, u, mat, p, cs, T);
  }

  //---------------------------------------------------------------------------
  // Define a new EOS here
  //---------------------------------------------------------------------------
  void idealGas(const ftype &rho, const ftype &rhom1, const ftype &u, 
		const int &mat, ftype &P, ftype &cs, ftype &T) const {
    ftype Rgas = 0.00820292998, gamma = 1.6666666, meanm=1.0;
    P  = (gamma - 1.0) * rho * u;
    cs = sqrt(gamma*P/rho);
    T  = meanm*(gamma-1.)*u/Rgas;
  }

  void tillotson(const ftype &rho, const ftype &rhom1, const ftype &u, 
		 const int &mat, ftype &P, ftype &cs, ftype &T) const {
    eosMat m = eosTab.get(mat);
    ftype PC     = 0.;
    ftype csC    = 0.;
    ftype rho0m1 = 1. / m.rho0;
    ftype eta    = rho*rho0m1;
    ftype mu     = eta - 1.;

    ftype csmin  = 0.25*m.A*rho0m1;
    ftype Pmin   = -1.e15;

    ftype c1     = u / (m.u0*eta*eta);
    ftype c2     = 1. / (c1 + 1.);

    if (u > m.Eiv && eta < 1.) {
      ftype d1  = m.rho0*rhom1;
      ftype d2  = d1 - 1.;
      ftype ex1 = exp(-m.beta*d2);
      ftype ex2 = exp(-m.alpha*d2*d2);

      P   = m.a*rho*u;
      P  += ex2*(m.b*rho*u*c2 + m.A*mu*ex1);
      
      cs  = m.b*u*(3.*c1 + 1)*c2*c2 + 2.*m.alpha*d2*m.b*d1*u*c2;
      cs += m.A*ex1*((2.*m.alpha*d2 + m.beta)*mu*d1*rhom1 + rho0m1);
      cs  = cs*ex2 + m.a*u;
      cs += P*rhom1*(m.a + m.b*c2*c2*ex2);
      if (cs < 0.) cs = 0.;
    }

    if (u < m.Ecv || eta >= 1.) {
      PC   = (m.a + m.b*c2)*rho*u + m.A*mu + m.B*mu*mu;

      csC  = m.a*u + rho0m1*(m.A + 2.*m.B*mu) + m.b*u*(3.*c1 + 1.)*c2*c2; 
      csC += PC*rhom1*(m.a + m.b*c2*c2);

      if (u > m.Eiv && u < m.Ecv && eta < 1) {
	ftype e1 = m.Ecv - u;
	ftype e2 = u - m.Eiv;
	ftype e3 = 1. / (m.Ecv - m.Eiv);
	P  = (e2 * P  + e1 * PC)  * e3;
  	cs = (e2 * cs + e1 * csC) * e3;
      } else {
	P  = PC;
  	cs = csC;
      }
    }
    if (cs < csmin) cs = csmin;
    if (P  < Pmin)  P  = Pmin;
    cs = sqrt(cs);
  }
  
#ifdef ANEOS
  void aneos(const ftype &rho, const ftype &rhom1, const ftype &u,
	     const int &mat, ftype &p, ftype &cs, ftype &t) const {
    int    KPA , MAT = mat;

    double RHO = (double)rho, U = (double)u, P = (double)p, CS = (double)cs,
           T   = (double)t;

    rooten(RHO, U, MAT, T, P, CS, KPA);
    p = (ftype)P; cs = (ftype)CS; t = (ftype)T;
  }

  void rooten(double &rhoi, double &ui, int &mati,
	      double &ti, double &pi, double &csi, int &kpai) const {
    const double eps   = 1.e-5;
    const int    itmax = 30;
    double       _S, _CV, _DPDT, _DPDR, _FKROS, _FME, _FMA;
    double       a, b, c, d, e, ei, fa, fb, fc, p, q, r, s, tm, tol1;

    // Initial temperature bracket (in eV)
    double t1 = 0.001, t2 = 6.0;

    // Check lower boundary
    for (fa = 0.0; fa >= 0.0; t1 *= 0.1) {
      if (t1 < 1.e-6) { ui = ei; ti = a; return; }
      a = t1;
      aneos_(&a, &rhoi, &pi, &ei, &_S, &_CV, &_DPDT, &_DPDR, &_FKROS, 
	     &csi, &kpai, &mati, &_FME, &_FMA);
      fa = ei - ui;      
    } 
      
    // Check upper boundary
    for (fb = 0.0; fb <= 0.0; t2 *= 3.0) {
      if (t2 > 1.e15) { std::cout << "Temperature out of bounds!"; exit(1); }
      b = t2;
      aneos_(&b, &rhoi, &pi, &ei, &_S, &_CV, &_DPDT, &_DPDR, &_FKROS, 
	     &csi, &kpai, &mati, &_FME, &_FMA);
      fb = ei - ui;      
    } 

    // Start iteration
    fc = fb; 
    
    // just to shut up compiler
    c = a; d = b - a; e = d; q = 0;
    
    for(int i = 0; i < itmax; i++) {
      if (fb*fc > 0)           { c = a; fc = fa; d = b - a; e = d; }
      if (fabs(fc) < fabs(fb)) { 
	a = b; b = c; c = a; fa = fb; fb = fc; fc = fa;
      }
      tm   = 0.5 * (c - b);
      tol1 = 2. * eps * fabs(b);
      if (fabs(tm) < tol1 || fabs(fb/ui) < eps) { ti = b; return; }
      if (fabs(e)  > tol1 && fabs(fa) > fabs(fb)) {
	s = fb / fa;
	if (a == c) { p = 2.*tm*s; } 
	else { 
	  q = fa/fc; r = fb/fc; 
	  p = s*(2.*tm*q*(q-r)-(b-a)*(r-1.));
	  q = (q-1.)*(r-1.)*(s-1.);
	}
	// there might be a problem with q here
	if (p > 0.) q = -q;
	p = fabs(p);
	if (2.*p < 3.*tm*q-fabs(tol1*q) && 2.*p < fabs(e*q)) { 
	  e = d; d = p/q; 
	}
	else { d = tm; e = d; }
      } else { d = tm; e = d; }
      a = b; fa = fb;
      if (fabs(d) > tol1) b += d; else {
	if (tm >= 0.) b += fabs(tol1); else b -= fabs(tol1);
      }
      aneos_(&b, &rhoi, &pi, &ei, &_S, &_CV, &_DPDT, &_DPDR, &_FKROS,
	     &csi, &kpai, &mati, &_FME, &_FMA);
      fb = ei - ui;
    }
    std::cout << "rooten: T iteration did not converge!" << std::endl;
    exit(1);
  }
#endif

  void findrho(const ftype &_T, const int &mat, ftype &rho0, ftype &_u0) const {
    switch (number) {
      // if your new EOS need to find density to minimize initial pressure,
      // add it here: case 3: ......
#ifdef ANEOS
    case 2: {
      double _S, _CV, _DPDT, _DPDR, _FKROS, _CS, _FME, _FMA;
      int    itmax = 50, _KPA, _MAT = 1;
      double a = rho0, b = rho0, c, dx = 0.1 * rho0, pa, pb, pc, tol = 1.e-5;
      double T = (double)_T, u0 = (double)_u0;
      do {
	a -= dx;
	aneos_(&T, &a, &pa, &u0, &_S, &_CV, &_DPDT, &_DPDR, &_FKROS,
	       &_CS, &_KPA, &_MAT, &_FME, &_FMA);
      } while (pa > 0. && a >= dx);
      do {
	b *= 1.5;
	aneos_(&T, &b, &pb, &u0, &_S, &_CV, &_DPDT, &_DPDR, &_FKROS, 
	       &_CS, &_KPA, &_MAT, &_FME, &_FMA);
      } while (pb < 0. && b < 30.*rho0);
      
      
      for (int i = 0; i < itmax && b - a > tol; i++) {
	c = a + (b - a) * .5;
	aneos_(&T, &c, &pc, &u0, &_S, &_CV, &_DPDT, &_DPDR, &_FKROS,
	       &_CS, &_KPA, &_MAT, &_FME, &_FMA);
	if (pc < 0.) { a = c; pa = pc; } else { b = c; pb = pc; }
      }
      rho0 = (ftype)c; _u0 = (ftype)u0;
      break;
    }
#endif
    default: { 
      int   itmax = 50;
      ftype a = rho0, b = rho0, c, dx = 0.001 * rho0;
      ftype pa, pb, pc, cs, tol = 1.e-10, T;
      do {
	a -= dx;
	eos(a, 1./a, _u0, mat, pa, cs, T);
      } while (pa > 0. && a > 0.9*rho0);
      do {
	b *= 1.005;
	eos(b, 1./b, _u0, mat, pb, cs, T);
      } while (pb < 0. && b < 0.1*rho0);
      
      for (int i = 0; i < itmax && b - a > tol; i++) {
	c = a + (b - a) * .5;
	eos(c, 1./c, _u0, mat, pc, cs, T);
	if (pc < 0.) { a = c; pa = pc; } else { b = c; pb = pc; }
      }
      rho0 = c;
      break;
    }
    }
  }
};

#endif
