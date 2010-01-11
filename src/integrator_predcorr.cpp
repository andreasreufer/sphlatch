#ifndef INTEGRATOR_PREDCORR_CPP
#define INTEGRATOR_PREDCORR_CPP

namespace sphlatch {
template<typename _T>
class PredictorCorrectorO1
{
public:

   void bootstrap(_T& _var, _T& _dvar)
   {
      odvar = _dvar;
   }

   void predict(_T& _var, _T& _dvar, const fType _dt)
   {
      ovar  = _var;
      _var += (1.5 * _dvar - 0.5 * odvar) * _dt;

      odvar = _dvar;
      //_dvar = 0.;
   }

   void correct(_T& _var, _T& _dvar, const fType _dt)
   {
      _var = ovar + 0.5 * _dt * (_dvar + odvar);
   }

   _T ovar, odvar;
};

template<typename _T>
class PredictorCorrectorO2
{
public:

   void bootstrap(_T& _var, _T& _dvar, _T& _ddvar)
   {
      odvar  = _dvar;
      oddvar = _ddvar;
   }

   void predict(_T& _var, _T& _dvar, _T& _ddvar, const fType _dt)
   {
      ovar  = _var;
      _var += (1.5 * _dvar - 0.5 * odvar) * _dt;

      odvar  = _dvar;
      _dvar += (1.5 * _ddvar - 0.5 * oddvar) * _dt;

      oddvar = _ddvar;
   }

   void correct(_T& _var, _T& _dvar, _T& _ddvar, const fType _dt)
   {
      _var  = ovar + 0.5 * _dt * (_dvar + odvar);
      _dvar = odvar + 0.5 * _dt * (_ddvar + oddvar);
   }

   _T ovar, odvar, oddvar;
};
};

#endif
