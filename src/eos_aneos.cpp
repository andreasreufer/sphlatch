#ifndef SPHLATCH_EOS_ANEOS
#define SPHLATCH_EOS_ANEOS

/*
 *  eos_aneos.h
 *
 *
 *  Created by Andreas Reufer on 15.09.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include <fstream>
#include <boost/lexical_cast.hpp>

#include "typedefs.h"
#include "eos_generic.cpp"
#include "err_handler.cpp"

#ifdef SPHLATCH_ANEOS_TABLE
 #include "lookup_table2D.cpp"
#endif

#ifdef SPHLATCH_HDF5
 #include "hdf5_io.cpp"
#endif

#include "brent.cpp"

///
/// FORTRAN subroutine "ANEOSINIT" to init ANEOS
///
extern "C" void aneosinit_(const char*,  /// materials file
                           int);         /// filename string length

///
/// FORTRAN subroutine "ANEOS" for p(rho,T,mat)
///
extern "C" void aneosd_(const double*,  /// T (input)
                        const double*,  /// rho (input)
                        double*,        /// p
                        double*,        /// E
                        double*,        /// S
                        double*,        /// c_v
                        double*,        /// dp/dT
                        double*,        /// dp/dRho
                        double*,        /// fkro (Rosseland mean opacity)
                        double*,        /// cs
                        int*,           /// phase
                        const int*,     /// material (input)
                        double*,        /// rhoL
                        double*         /// rhoH
                        );

extern "C" void aneosv_(const int*,     /// no of array elements
                        const double*,  /// T (input)
                        const double*,  /// rho (input)
                        const int*,     /// material (input)
                        double*,        /// p
                        double*,        /// E
                        double*,        /// S
                        double*,        /// c_v
                        double*,        /// dp/dT
                        double*,        /// dp/dRho
                        double*,        /// fkro (Rosseland mean opacity)
                        double*,        /// cs
                        int*,           /// phase
                        double*,        /// lower  phase density
                        double*,        /// higher phase density
                        double*         /// ionization number
                        );

namespace sphlatch {
template<typename lutT>
class ANEOStables {
   // FIXME: get ranges
public:
#ifdef SPHLATCH_HDF5
   ANEOStables<lutT>(const std::string _fname, const std::string _path,
                     const std::string _xname, const std::string _yname,
                     const std::string _zname) :
      T(_fname, _path + _xname, _path + _yname, _path + "/T"),
      p(_fname, _path + _xname, _path + _yname, _path + "/p"),
      cv(_fname, _path + _xname, _path + _yname, _path + "/cv"),
      dpdt(_fname, _path + _xname, _path + _yname, _path + "/dpdt"),
      dpdrho(_fname, _path + _xname, _path + _yname, _path + "/dpdrho"),
      fkros(_fname, _path + _xname, _path + _yname, _path + "/fkros"),
      cs(_fname, _path + _xname, _path + _yname, _path + "/cs"),
      rhoL(_fname, _path + _xname, _path + _yname, _path + "/rhoL"),
      rhoH(_fname, _path + _xname, _path + _yname, _path + "/rhoH"),
      phase(_fname, _path + _xname, _path + _yname, _path + "/phase"),
      z(_fname, _path + _xname, _path + _yname, _path + _zname),
      xname(_xname),
      yname(_yname),
      zname(_zname)
   { };
#endif
   ANEOStables<lutT>(const fvectT &_logx, const fvectT &_logy,
                     const fmatrT &_T,
                     const fmatrT &_p,
                     const fmatrT &_cv,
                     const fmatrT &_dpdt,
                     const fmatrT &_dpdrho,
                     const fmatrT &_fkros,
                     const fmatrT &_cs,
                     const fmatrT &_rhoL,
                     const fmatrT &_rhoH,
                     const fmatrT &_phase,
                     const fmatrT &_z,
                     const std::string _xname,
                     const std::string _yname,
                     const std::string _zname) :
      T(_logx, _logy, _T),
      p(_logx, _logy, _p),
      cv(_logx, _logy, _cv),
      dpdt(_logx, _logy, _dpdt),
      dpdrho(_logx, _logy, _dpdrho),
      fkros(_logx, _logy, _fkros),
      cs(_logx, _logy, _cs),
      rhoL(_logx, _logy, _rhoL),
      rhoH(_logx, _logy, _rhoH),
      phase(_logx, _logy, _phase),
      z(_logx, _logy, _z),
      xname(_xname),
      yname(_yname),
      zname(_zname)
   { };

   ~ANEOStables<lutT>() { };

#ifdef SPHLATCH_HDF5
   void storeTable(const std::string _fname, const std::string _path)
   {
      HDF5File tf(_fname);

      tf.doublePrecOut();

      tf.createGroup(_path);
      tf.savePrimitive(_path + xname, T.getX());
      tf.savePrimitive(_path + yname, T.getY());
      tf.savePrimitive(_path + zname, z.getF());

      tf.savePrimitive(_path + "/T", T.getF());
      tf.savePrimitive(_path + "/p", p.getF());
      tf.savePrimitive(_path + "/cv", cv.getF());
      tf.savePrimitive(_path + "/dpdt", dpdt.getF());
      tf.savePrimitive(_path + "/dpdrho", dpdrho.getF());
      tf.savePrimitive(_path + "/fkros", fkros.getF());
      tf.savePrimitive(_path + "/cs", cs.getF());
      tf.savePrimitive(_path + "/rhoL", rhoL.getF());
      tf.savePrimitive(_path + "/rhoH", rhoH.getF());
      tf.savePrimitive(_path + "/phase", phase.getF());
   }
#endif


private:
   lutT        T, p, cv, dpdt, dpdrho, fkros, cs, rhoL, rhoH, phase, z;
   std::string xname, yname, zname;

public:
   void operator()(const fType _x, const fType _y,
                   fType& _T,
                   fType& _p,
                   fType& _cv,
                   fType& _dpdt,
                   fType& _dpdrho,
                   fType& _fkros,
                   fType& _cs,
                   fType& _rhoL,
                   fType& _rhoH,
                   iType& _phase,
                   fType& _z)
   {
      const fType lx = log10(_x);
      const fType ly = log10(_y);

      _p = p(lx, ly);
      const size_t ixl = p.ixl;
      const size_t iyl = p.iyl;

      _T      = T(lx, ly, ixl, iyl);
      _p      = p(lx, ly, ixl, iyl);
      _cv     = cv(lx, ly, ixl, iyl);
      _dpdt   = dpdt(lx, ly, ixl, iyl);
      _dpdrho = dpdrho(lx, ly, ixl, iyl);
      _fkros  = fkros(lx, ly, ixl, iyl);
      _cs     = cs(lx, ly, ixl, iyl);
      _rhoL   = rhoL(lx, ly, ixl, iyl);
      _rhoH   = rhoH(lx, ly, ixl, iyl);
      _phase  = lrint(phase(lx, ly, ixl, iyl));
      _z      = z(lx, ly, ixl, iyl);
   }

   // just forward the getRange() method
   void getRange(fType& _xMin, fType& _xMax, fType& _yMin, fType& _yMax)
   {
      p.getRange(_xMin, _xMax, _yMin, _yMax);
   }
};



template<typename _partT>
class ANEOS : public EOS {
public:
   ANEOS();
   ~ANEOS();

   static ANEOS& instance();

   static ANEOS* _instance;

public:
   void operator()(_partT& _part);

   void rootS(fType& _T, const fType _rho, const iType _mat,
              fType& _p, fType& _u, const fType _S, fType& _cv,
              fType& _dpdt, fType& _dpdrho, fType& _fkros, fType& _cs,
              iType& _kpa, fType& _rhoL, fType& _rhoH);

   void rootU(fType& _T, const fType _rho, const iType _mat,
              fType& _p, const fType _u, fType& _S, fType& _cv,
              fType& _dpdt, fType& _dpdrho, fType& _fkros, fType& _cs,
              iType& _kpa, fType& _rhoL, fType& _rhoH);

   void rootP(const fType _T, fType& _rho, const iType _mat,
              const fType _p, fType& _u, fType& _S, fType& _cv,
              fType& _dpdt, fType& _dpdrho, fType& _fkros, fType& _cs,
              iType& _kpa, fType& _rhoL, fType& _rhoH);

   static void callaneos(const fType _T, const fType _rho, const iType _mat,
                         fType& _p, fType& _u, fType& _S, fType& _cv,
                         fType& _dpdt, fType& _dpdrho, fType& _fkros,
                         fType& _cs,
                         identType& _kpa, fType& _rhoL,
                         fType& _rhoH);

#ifdef SPHLATCH_ANEOS_TABLE
   void tableU(fType& _T, const fType _rho, const iType _mat,
               fType& _p, const fType _u, fType& _S, fType& _cv,
               fType& _dpdt, fType& _dpdrho, fType& _fkros, fType& _cs,
               iType& _kpa, fType& _rhoL, fType& _rhoH);

   void tableS(fType& _T, const fType _rho, const iType _mat,
               fType& _p, fType& _u, const fType _S, fType& _cv,
               fType& _dpdt, fType& _dpdrho, fType& _fkros, fType& _cs,
               iType& _kpa, fType& _rhoL, fType& _rhoH);

 #ifdef SPHLATCH_HDF5
   void loadTableU(std::string _fname, const identType _mat);
   void loadTableS(std::string _fname, const identType _mat);

   void storeTableU(std::string _fname, const identType _mat);
   void storeTableS(std::string _fname, const identType _mat);

private:
   std::string getHDFPath(const iType _mat);
 #endif
#endif




private:
   fType nan;
   bool  quiet;

   struct uRoot;
   struct SRoot;
   struct pRoot;

   typedef BrentRooter<uRoot>   uRooterT;
   typedef BrentRooter<SRoot>   SRooterT;
   typedef BrentRooter<pRoot>   pRooterT;

   uRooterT uRooter;
   SRooterT SRooter;
   pRooterT pRooter;

#ifdef SPHLATCH_ANEOS_TABLE
   const size_t        maxMatId;
   fvectT              rhoMin, rhoMed, rhoMax, uMin, uMax, SMin, SMax;
   std::vector<size_t> npRhoL, npRhoH, npU, npS;

   typedef LookupTable2D<InterpolateBilinear>   lutT;
   std::vector<ANEOStables<lutT> *> tblUL, tblUH, tblSL, tblSH;
   std::vector<bool>                tblUinit, tblSinit;

   template<typename lutT, typename rooterT>
   ANEOStables<lutT> * generateTable(const iType _mat, const size_t _npx,
                                     const size_t _npy, const fType _xmin,
                                     const fType _xmax, const fType _ymin,
                                     const fType _ymax, rooterT& _rooter,
                                     const std::string _xname,
                                     const std::string _yname,
                                     const std::string _zname);
#endif

};


template<typename _partT>
ANEOS<_partT>::ANEOS()
#ifdef SPHLATCH_ANEOS_TABLE
   : maxMatId(32),
     rhoMin(maxMatId + 1),
     rhoMed(maxMatId + 1),
     rhoMax(maxMatId + 1),
     uMin(maxMatId + 1),
     uMax(maxMatId + 1),
     SMin(maxMatId + 1),
     SMax(maxMatId + 1),
     npRhoL(maxMatId + 1),
     npRhoH(maxMatId + 1),
     npU(maxMatId + 1),
     npS(maxMatId + 1),
     tblUL(maxMatId + 1),
     tblUH(maxMatId + 1),
     tblSL(maxMatId + 1),
     tblSH(maxMatId + 1),
     tblUinit(maxMatId + 1),
     tblSinit(maxMatId + 1)
#endif
{
#ifdef SPHLATCH_MANEOS
   std::string matFilename = "ANEOS.QUARTZ";
#else
   std::string matFilename = "aneos.input";
#endif

   aneosinit_(matFilename.c_str(), matFilename.size());

#ifdef SPHLATCH_LOGGER
   Logger.stream << "init ANEOS EOS with file "
                 << matFilename;
   Logger.flushStream();
#endif

#ifdef SPHLATCH_ANEOS_TABLE
   for (size_t i = 0; i < maxMatId + 1; i++)
   {
      npRhoL[i] = 100;
      npRhoH[i] = 500;

      npU[i] = 1500;
      npS[i] = 1500;

      uMin[i] = 1.e7;   // erg/g
      uMax[i] = 1.e14;

      SMin[i] = 1.e9;   // erg/eV/g
      SMax[i] = 2.e13;

      rhoMin[i] = 1.0e-10; // g/cm3
      rhoMed[i] = 0.1;
      rhoMax[i] = 25.;

      tblUL[i] = NULL;
      tblUH[i] = NULL;

      tblSL[i] = NULL;
      tblSH[i] = NULL;

      tblUinit[i] = false;
      tblSinit[i] = false;
   }

   // iron mat 5
   uMin[5] = 1.e8;
   uMax[5] = 1.e14;

   // H2O mat 2
   rhoMax[2] = 7.; // g/cm3

   // SiO2 mat 1
   uMin[1]   = 1.e8;
   uMax[1]   = 1.e16;
   rhoMax[1] = 50.;
#endif


   nan = std::numeric_limits<fType>::quiet_NaN();

   ///
   /// throw exceptions per default, keep quiet when generating tables
   ///
   quiet = false;
}

template<typename _partT>
ANEOS<_partT>::~ANEOS() { }

template<typename _partT>
ANEOS<_partT> * ANEOS<_partT>::_instance = NULL;

template<typename _partT>
ANEOS<_partT>& ANEOS<_partT>::instance()
{
   if (_instance == NULL)
      _instance = new ANEOS;
   return(*_instance);
}

///
/// common EOS interface for particle use
///
template<typename _partT>
void ANEOS<_partT>::operator()(_partT& _part)
{
   fType dummy;

#ifdef SPHLATCH_ANEOS_TABLE
   const iType mat = _part.mat;
   // fall back to rooting algorithm, if outside the table values
   if (
 #ifdef SPHLATCH_TIMEDEP_ENTROPY
      (_part.S < SMin[mat]) or (_part.S > SMax[mat])
 #else
      (_part.u < uMin[mat]) or (_part.u > uMax[mat])
 #endif
      or (_part.rho < rhoMin[mat]) or (_part.rho > rhoMax[mat]))
   {
#endif

#ifdef SPHLATCH_TIMEDEP_ENTROPY
   rootS(_part.T, _part.rho, _part.mat, _part.p, _part.u, _part.S, dummy,
         dummy, dummy, dummy, _part.cs, _part.phase, _part.rhoL, _part.rhoH);
#else
   rootU(_part.T, _part.rho, _part.mat, _part.p, _part.u, _part.S, dummy,
         dummy, dummy, dummy, _part.cs, _part.phase, _part.rhoL, _part.rhoH);
#endif

#ifdef SPHLATCH_ANEOS_TABLE
}

else
{
 #ifdef SPHLATCH_TIMEDEP_ENTROPY
   tableS(_part.T, _part.rho, _part.mat, _part.p, _part.u, _part.S, dummy,
          dummy, dummy, dummy, _part.cs, _part.phase, _part.rhoL, _part.rhoH);
 #else
   tableU(_part.T, _part.rho, _part.mat, _part.p, _part.u, _part.S, dummy,
          dummy, dummy, dummy, _part.cs, _part.phase, _part.rhoL, _part.rhoH);
 #endif
}
#endif


#ifdef SPHLATCH_NONEGPRESS
   if (_part.p < 0.)
      _part.p = 0.;
#endif
   if (_part.p != _part.p)
   {
#ifdef SPHLATCH_LOGGER
      Logger.stream << "part ID: " << _part.id << " rho: " << _part.rho <<
      " u: " << _part.u << " mat: " << _part.mat << " has NaN pressure!";
      Logger.flushStream();
#endif
   }
}

template<typename _partT>
void ANEOS<_partT>::rootS(fType& _T, const fType _rho, const iType _mat,
                          fType& _p, fType& _u, const fType _S, fType& _cv,
                          fType& _dpdt, fType& _dpdrho, fType& _fkros,
                          fType& _cs,
                          iType& _kpa, fType& _rhoL,
                          fType& _rhoH)
{
   fType Tmin = 0.001;
   fType Tmax = 6.;

   SRooter.f.x     = _rho;
   SRooter.f.mat   = _mat;
   SRooter.f.ytarg = _S;

   while (SRooter.f(Tmin) > 0.)
      Tmin *= 0.1;

   while (SRooter.f(Tmax) < 0.)
      Tmax *= 3.0;

   _T      = SRooter(Tmin, Tmax, 1.e-5);
   _p      = SRooter.f.p;
   _u      = SRooter.f.z;
   _cv     = SRooter.f.cv;
   _dpdt   = SRooter.f.dpdt;
   _dpdrho = SRooter.f.dpdrho;
   _fkros  = SRooter.f.fkros;
   _cs     = SRooter.f.cs;
   _kpa    = SRooter.f.kpa;
   _rhoL   = SRooter.f.rhoL;
   _rhoH   = SRooter.f.rhoH;
}

template<typename _partT>
void ANEOS<_partT>::rootU(fType& _T, const fType _rho, const iType _mat,
                          fType& _p, const fType _u, fType& _S, fType& _cv,
                          fType& _dpdt, fType& _dpdrho, fType& _fkros,
                          fType& _cs,
                          iType& _kpa, fType& _rhoL,
                          fType& _rhoH)
{
   fType Tmin = 0.001;
   fType Tmax = 6.;

   uRooter.f.x     = _rho;
   uRooter.f.mat   = _mat;
   uRooter.f.ytarg = _u;

   while (uRooter.f(Tmin) > 0.)
      Tmin *= 0.1;

   while (uRooter.f(Tmax) < 0.)
      Tmax *= 3.0;

   //std::cout << _rho << " " << _S << " " << Tmin << " " << Tmax << "\n";
   _T      = uRooter(Tmin, Tmax, 1.e-5);
   _p      = uRooter.f.p;
   _S      = uRooter.f.z;
   _cv     = uRooter.f.cv;
   _dpdt   = uRooter.f.dpdt;
   _dpdrho = uRooter.f.dpdrho;
   _fkros  = uRooter.f.fkros;
   _cs     = uRooter.f.cs;
   _kpa    = uRooter.f.kpa;
   _rhoL   = uRooter.f.rhoL;
   _rhoH   = uRooter.f.rhoH;
}

template<typename _partT>
void ANEOS<_partT>::rootP(const fType _T, fType& _rho, const iType _mat,
                          const fType _p, fType& _u, fType& _S, fType& _cv,
                          fType& _dpdt, fType& _dpdrho, fType& _fkros,
                          fType& _cs,
                          iType& _kpa, fType& _rhoL,
                          fType& _rhoH)
{
   fType rhoMin = 1.e-3;
   fType rhoMax = 10.;

   pRooter.f.T     = _T;
   pRooter.f.mat   = _mat;
   pRooter.f.ptarg = _p;

   while (pRooter.f(rhoMin) > 0.)
      rhoMin *= 0.1;

   while (pRooter.f(rhoMax) < 0.)
      rhoMax *= 3.0;

   _rho    = pRooter(rhoMin, rhoMax, 1.e-5);
   _u      = pRooter.f.u;
   _S      = pRooter.f.S;
   _cv     = pRooter.f.cv;
   _dpdt   = pRooter.f.dpdt;
   _dpdrho = pRooter.f.dpdrho;
   _fkros  = pRooter.f.fkros;
   _cs     = pRooter.f.cs;
   _kpa    = pRooter.f.kpa;
   _rhoL   = pRooter.f.rhoL;
   _rhoH   = pRooter.f.rhoH;
}

template<typename _partT>
void ANEOS<_partT>::callaneos(const fType _T, const fType _rho,
                              const iType _mat,
                              fType& _p, fType& _u, fType& _S, fType& _cv,
                              fType& _dpdt, fType& _dpdrho, fType& _fkros,
                              fType& _cs,
                              identType& _kpa, fType& _rhoL,
                              fType& _rhoH)
{
   static double T, rho, p, u, S, cv, dpdt, dpdrho, fkros, cs, rhoL, rhoH, ion;
   static int    kpa, mat;

   const int n = 1;

   rho = static_cast<double>(_rho);
   T   = static_cast<double>(_T);
   mat = static_cast<int>(_mat);

   aneosv_(&n, &T, &rho, &mat, &p, &u, &S, &cv, &dpdt, &dpdrho, &fkros,
           &cs, &kpa, &rhoL, &rhoH, &ion);

   _p      = static_cast<fType>(p);
   _u      = static_cast<fType>(u);
   _S      = static_cast<fType>(S);
   _cv     = static_cast<fType>(cv);
   _dpdt   = static_cast<fType>(dpdt);
   _dpdrho = static_cast<fType>(dpdrho);
   _fkros  = static_cast<fType>(fkros);
   _cs     = static_cast<fType>(cs);
   _kpa    = static_cast<identType>(kpa);
   _rhoL   = static_cast<fType>(rhoL);
   _rhoH   = static_cast<fType>(rhoH);
}

#ifdef SPHLATCH_ANEOS_TABLE
template<typename _partT>
void ANEOS<_partT>::tableU(fType& _T, const fType _rho, const iType _mat,
                           fType& _p, const fType _u, fType& _S, fType& _cv,
                           fType& _dpdt, fType& _dpdrho, fType& _fkros,
                           fType& _cs,
                           iType& _kpa, fType& _rhoL, fType& _rhoH)
{
   if (not tblUinit[_mat])
   {
      tblUL[_mat] = generateTable<lutT, uRooterT>(_mat, npRhoL[_mat], npU[_mat],
                                                  rhoMin[_mat], rhoMed[_mat],
                                                  uMin[_mat], uMax[_mat],
                                                  uRooter,
                                                  "/log10rho", "/log10u",
                                                  "/S");
      tblUH[_mat] = generateTable<lutT, uRooterT>(_mat, npRhoH[_mat], npU[_mat],
                                                  rhoMed[_mat], rhoMax[_mat],
                                                  uMin[_mat], uMax[_mat],
                                                  uRooter,
                                                  "/log10rho", "/log10u",
                                                  "/S");
      tblUinit[_mat] = true;
   }

   if (_rho < rhoMed[_mat])
      tblUL[_mat]->operator()(_rho, _u, _T, _p, _cv, _dpdt, _dpdrho, _fkros,
                              _cs, _rhoL, _rhoH, _kpa, _S);
   else
      tblUH[_mat]->operator()(_rho, _u, _T, _p, _cv, _dpdt, _dpdrho, _fkros,
                              _cs, _rhoL, _rhoH, _kpa, _S);
}

template<typename _partT>
void ANEOS<_partT>::tableS(fType& _T, const fType _rho, const iType _mat,
                           fType& _p, fType& _u, const fType _S, fType& _cv,
                           fType& _dpdt, fType& _dpdrho, fType& _fkros,
                           fType& _cs,
                           iType& _kpa, fType& _rhoL, fType& _rhoH)
{
   if (not tblSinit[_mat])
   {
      tblSL[_mat] = generateTable<lutT, SRooterT>(_mat, npRhoL[_mat], npS[_mat],
                                                  rhoMin[_mat], rhoMed[_mat],
                                                  SMin[_mat], SMax[_mat],
                                                  SRooter,
                                                  "/log10rho", "/log10S",
                                                  "/u");
      tblSH[_mat] = generateTable<lutT, SRooterT>(_mat, npRhoH[_mat], npS[_mat],
                                                  rhoMed[_mat], rhoMax[_mat],
                                                  SMin[_mat], SMax[_mat],
                                                  SRooter,
                                                  "/log10rho", "/log10S",
                                                  "/u");
      tblSinit[_mat] = true;
   }

   if (_rho < rhoMed[_mat])
      tblSL[_mat]->operator()(_rho, _S, _T, _p, _cv, _dpdt, _dpdrho, _fkros,
                              _cs, _rhoL, _rhoH, _kpa, _u);
   else
      tblSH[_mat]->operator()(_rho, _S, _T, _p, _cv, _dpdt, _dpdrho, _fkros,
                              _cs, _rhoL, _rhoH, _kpa, _u);
}

 #ifdef SPHLATCH_HDF5
template<typename _partT>
void ANEOS<_partT>::loadTableU(std::string _fname, const identType _mat)
{
   if (tblUinit[_mat])
      delete tblUL[_mat], tblUH[_mat];

   std::string basepath = getHDFPath(_mat);
   tblUL[_mat] =
      new ANEOStables<lutT>(_fname, basepath + "UL", "/log10rho", "/log10u",
                            "/S");


   tblUH[_mat] =
      new ANEOStables<lutT>(_fname, basepath + "UH", "/log10rho", "/log10u",
                            "/S");

   fType logRhoMin, logRhoMed, logRhoMax, logUMin, logUMax;
   tblUL[_mat]->getRange(logRhoMin, logRhoMed, logUMin, logUMax);
   tblUH[_mat]->getRange(logRhoMed, logRhoMax, logUMin, logUMax);

   // override limits
   rhoMin[_mat] = pow(10., logRhoMin);
   rhoMed[_mat] = pow(10., logRhoMed);
   rhoMax[_mat] = pow(10., logRhoMax);

   uMin[_mat] = pow(10., logUMin);
   uMax[_mat] = pow(10., logUMax);

  #ifdef SPHLATCH_LOGGER
   Logger.stream << "ANEOS U tables loaded for mat id " << _mat << "\n"
                 << "   rho: " << rhoMin[_mat]
                 << " - " << rhoMed[_mat]
                 << " - " << rhoMax[_mat] << "\n"
                 << "   u:   " << uMin[_mat]
                 << " - " << uMax[_mat] << "\n";
   Logger.flushStream();
  #endif

   tblUinit[_mat] = true;
}

template<typename _partT>
void ANEOS<_partT>::loadTableS(std::string _fname, const identType _mat)
{
   if (tblSinit[_mat])
      delete tblSL[_mat], tblSH[_mat];

   std::string basepath = getHDFPath(_mat);
   tblSL[_mat] =
      new ANEOStables<lutT>(_fname, basepath + "SL", "/log10rho", "/log10S",
                            "/u");
   tblSH[_mat] =
      new ANEOStables<lutT>(_fname, basepath + "SH", "/log10rho", "/log10S",
                            "/u");

   fType logRhoMin, logRhoMed, logRhoMax, logSMin, logSMax;
   tblSL[_mat]->getRange(logRhoMin, logRhoMed, logSMin, logSMax);
   tblSH[_mat]->getRange(logRhoMed, logRhoMax, logSMin, logSMax);

   // override limits
   rhoMin[_mat] = pow(10., logRhoMin);
   rhoMed[_mat] = pow(10., logRhoMed);
   rhoMax[_mat] = pow(10., logRhoMax);

   SMin[_mat] = pow(10., logSMin);
   SMax[_mat] = pow(10., logSMax);


  #ifdef SPHLATCH_LOGGER
   Logger.stream << "ANEOS S tables loaded for mat id " << _mat << "\n"
                 << "   rho: " << rhoMin[_mat]
                 << " - " << rhoMed[_mat]
                 << " - " << rhoMax[_mat] << "\n"
                 << "   S:   " << SMin[_mat]
                 << " - " << SMax[_mat] << "\n";
   Logger.flushStream();
  #endif

   tblSinit[_mat] = true;
}

template<typename _partT>
void ANEOS<_partT>::storeTableU(std::string _fname, const identType _mat)
{
   if (not tblUinit[_mat])
   {
      tblUL[_mat] = generateTable<lutT, uRooterT>(_mat, npRhoL[_mat], npU[_mat],
                                                  rhoMin[_mat], rhoMed[_mat],
                                                  uMin[_mat], uMax[_mat],
                                                  uRooter,
                                                  "/log10rho", "/log10u",
                                                  "/S");
      tblUH[_mat] = generateTable<lutT, uRooterT>(_mat, npRhoH[_mat], npU[_mat],
                                                  rhoMed[_mat], rhoMax[_mat],
                                                  uMin[_mat], uMax[_mat],
                                                  uRooter,
                                                  "/log10rho", "/log10u",
                                                  "/S");
      tblUinit[_mat] = true;
   }

   std::string basepath = getHDFPath(_mat);
   tblUL[_mat]->storeTable(_fname, basepath + "UL");
   tblUH[_mat]->storeTable(_fname, basepath + "UH");
}

template<typename _partT>
void ANEOS<_partT>::storeTableS(std::string _fname, const identType _mat)
{
   if (not tblSinit[_mat])
   {
      tblSL[_mat] = generateTable<lutT, SRooterT>(_mat, npRhoL[_mat], npS[_mat],
                                                  rhoMin[_mat], rhoMed[_mat],
                                                  SMin[_mat], SMax[_mat],
                                                  SRooter,
                                                  "/log10rho", "/log10S",
                                                  "/u");
      tblSH[_mat] = generateTable<lutT, SRooterT>(_mat, npRhoH[_mat], npS[_mat],
                                                  rhoMed[_mat], rhoMax[_mat],
                                                  SMin[_mat], SMax[_mat],
                                                  SRooter,
                                                  "/log10rho", "/log10S",
                                                  "/u");
      tblSinit[_mat] = true;
   }

   std::string basepath = getHDFPath(_mat);
   tblSL[_mat]->storeTable(_fname, basepath + "SL");
   tblSH[_mat]->storeTable(_fname, basepath + "SH");
}

template<typename _partT>
std::string ANEOS<_partT>::getHDFPath(const identType _mat)
{
   std::stringstream matstr, basepath;

   matstr << _mat;
   basepath << "/mat";
   for (size_t i = matstr.str().size(); i < 2; i++)
      basepath << "0";
   basepath << _mat;

   return(basepath.str());
}
 #endif




template<typename _partT>
template<typename lutT, typename rooterT>
ANEOStables<lutT> * ANEOS<_partT>::generateTable(const iType       _mat,
                                                 const size_t      _npx,
                                                 const size_t      _npy,
                                                 const fType       _xmin,
                                                 const fType       _xmax,
                                                 const fType       _ymin,
                                                 const fType       _ymax,
                                                 rooterT&          _rooter,
                                                 const std::string _xname,
                                                 const std::string _yname,
                                                 const std::string _zname)
{
 #ifdef SPHLATCH_LOGGER
   Logger.stream << "ANEOS generate table for mat id " << _mat;
   Logger.flushStream();
 #endif
   fvectT      logx(_npx);
   const fType dlogx   = (log10(_xmax) - log10(_xmin)) / (_npx - 1);
   const fType logxmin = log10(_xmin);

   for (size_t i = 0; i < _npx; i++)
      logx(i) = logxmin + dlogx * static_cast<fType>(i);

   fvectT      logy(_npy);
   const fType dlogy   = (log10(_ymax) - log10(_ymin)) / (_npy - 1);
   const fType logymin = log10(_ymin);

   for (size_t j = 0; j < _npy; j++)
      logy(j) = logymin + dlogy * static_cast<fType>(j);

   fmatrT T(_npx, _npy);
   fmatrT p(_npx, _npy);
   fmatrT cv(_npx, _npy);
   fmatrT dpdt(_npx, _npy);
   fmatrT dpdrho(_npx, _npy);
   fmatrT fkros(_npx, _npy);
   fmatrT cs(_npx, _npy);
   fmatrT phase(_npx, _npy);
   fmatrT z(_npx, _npy);
   fmatrT rhoL(_npx, _npy);
   fmatrT rhoH(_npx, _npy);

   const fType eps = 1.e-5;

   _rooter.quiet = true;
   _rooter.f.mat = _mat;

   size_t badpts = 0;

   for (size_t i = 0; i < _npx; i++)
   {
      _rooter.f.x = pow(10., logx(i));

      for (size_t j = 0; j < _npy; j++)
      {
         _rooter.f.ytarg = pow(10., logy(j));

         // use fixed Tmin and Tmax values for iron, comment the next two loops
         //const fType Tmin = 1.e-6;
         //const fType Tmax = 100.0;

         fType Tmin = 1.e-5;
         while (_rooter.f(Tmin) > 0.)
            Tmin *= 0.1;

         fType Tmax = 5.;
         while (_rooter.f(Tmax) < 0.)
            Tmax *= 3.0;

         const fType Ti = _rooter(Tmin, Tmax, eps);
         if (Ti == Ti)
         {
            T(i, j)      = Ti;
            p(i, j)      = _rooter.f.p;
            cv(i, j)     = _rooter.f.cv;
            dpdt(i, j)   = _rooter.f.dpdt;
            dpdrho(i, j) = _rooter.f.dpdrho;
            fkros(i, j)  = _rooter.f.fkros;
            cs(i, j)     = _rooter.f.cs;
            z(i, j)      = _rooter.f.z;
            rhoL(i, j)   = _rooter.f.rhoL;
            rhoH(i, j)   = _rooter.f.rhoH;
            phase(i, j)  = _rooter.f.kpa;
         }
         else
         {
            T(i, j)      = fTypeNan;
            p(i, j)      = fTypeNan;
            cv(i, j)     = fTypeNan;
            dpdt(i, j)   = fTypeNan;
            dpdrho(i, j) = fTypeNan;
            fkros(i, j)  = fTypeNan;
            cs(i, j)     = fTypeNan;
            z(i, j)      = fTypeNan;
            rhoL(i, j)   = fTypeNan;
            rhoH(i, j)   = fTypeNan;
            phase(i, j)  = fTypeNan;
            badpts++;
         }
      }
   }
   _rooter.quiet = false;

 #ifdef SPHLATCH_LOGGER
   Logger.stream << "ANEOS table (" << _npx << "x" << _npy
                 << ") initialised for mat id " << _mat
                 << " (range " << _xname << ": " << _xmin << "..." << _xmax
                 << ", range " << _yname << ": " << _ymin << "..." << _ymax
                 << ", " << badpts
                 << "/" << _npx * _npy
                 << " bad points)";
   Logger.flushStream();
 #endif

   ANEOStables<lutT> * tblPtr;
   tblPtr = new ANEOStables<lutT>(logx, logy, T, p, cv, dpdt, dpdrho, fkros,
                                  cs, rhoL, rhoH, phase, z, _xname, _yname,
                                  _zname);
   return(tblPtr);
}
#endif




template<typename _partT>
struct ANEOS<_partT>::uRoot
{
   fType operator()(const fType _T)
   {
      T = _T;

      callaneos(T, x, mat, p, y, z, cv, dpdt, dpdrho, fkros, cs, kpa, rhoL,
                rhoH);
      return(y - ytarg);
   }

   fType T, p, cv, dpdt, dpdrho, fkros, cs, rhoL, rhoH;
   iType mat, kpa;
   fType x, y, z, ytarg;
};

template<typename _partT>
struct ANEOS<_partT>::SRoot
{
   fType operator()(const fType _T)
   {
      T = _T;

      callaneos(T, x, mat, p, z, y, cv, dpdt, dpdrho, fkros, cs, kpa, rhoL,
                rhoH);
      return(y - ytarg);
   }

   fType T, p, cv, dpdt, dpdrho, fkros, cs, rhoL, rhoH;
   iType mat, kpa;
   fType Starg;
   fType x, y, z, ytarg;
};

template<typename _partT>
struct ANEOS<_partT>::pRoot
{
   // f(x) to find the root
   fType operator()(const fType _rho)
   {
      rho = _rho;

      callaneos(T, rho, mat, p, u, S, cv, dpdt, dpdrho, fkros, cs, kpa, rhoL,
                rhoH);
      return(p - ptarg);
   }

   fType T, rho, p, u, S, cv, dpdt, dpdrho, fkros, cs, rhoL, rhoH;
   iType mat, kpa;
   fType ptarg;
};
}
#endif
