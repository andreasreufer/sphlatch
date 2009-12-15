#ifndef SPHLATCH_IMPACT_SOLVER_CPP
#define SPHLATCH_IMPACT_SOLVER_CPP

/*
 *  impact_solver.cpp
 *
 *
 *  Created by Andreas Reufer on 10.12.09
 *  Copyright 2009 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"

namespace sphlatch {
void getImpactSetup(
   const fType _m_tar,
   const fType _m_imp,
   const fType _R_tar,
   const fType _R_imp,
   const fType _v_inf,
   const fType _L_tot,
   const fType _r_0,
   const fType _G,
   vect3dT&    _r_tarv,
   vect3dT&    _r_impv,
   vect3dT&    _v_tarv,
   vect3dT&    _v_impv,
   fType&      _t_0,
   fType&      _b_scal
   )
{
   const fType pi = M_PI;
   const fType m_tot = _m_tar + _m_imp;
   const fType gamma = _m_imp / m_tot;

   const fType v_esc = sqrt(2. * _G * m_tot / (_R_tar + _R_imp));
   const fType v_imp = sqrt(v_esc * v_esc + _v_inf * _v_inf);
   const fType r_imp = _R_tar + _R_imp;

   const fType mu = _G * m_tot;
   const fType k1 = r_imp * (v_imp * v_imp) / mu;

   const fType v_graz = sqrt((2. * mu / r_imp) + _v_inf * _v_inf);
   const fType L_graz = v_graz * _m_imp * (_R_tar + _R_imp);

   _b_scal   = _L_tot / L_graz;
   const fType beta_imp = acos(_b_scal);

   const fType E_tot = (v_imp * v_imp / 2.) - (mu / r_imp);

   const fType e = sqrt(pow((k1 - 1.) * cos(beta_imp), 2.)
                        + pow(sin(beta_imp), 2.));

   const fType theta_imp = atan2(k1 * sin(beta_imp) * cos(beta_imp),
                                 k1 * cos(beta_imp) * cos(beta_imp) - 1.);

   const fType r_per = r_imp * (1. + e * cos(theta_imp)) / (1. + e);
   const fType v_per = sqrt(pow(_v_inf, 2.) + (2. * mu / r_per));

   const fType v_0     = sqrt(_v_inf * _v_inf + (2. * mu / _r_0));
   const fType theta_0 = acos((r_per * (1. + e) - _r_0) / (e * _r_0));
   const fType beta_0  = acos(_L_tot / _r_0 * v_0 * _m_imp);

   fType k2, t_imp;

   if (e > 1.)
   {
      const fType a = r_per / (e - 1.);
      k2 = sqrt(mu / pow(a, 3.));

      t_imp = ((e * sqrt(e * e - 1.) * sin(theta_imp)) /
               (1. + e * cos(theta_imp)) -
               log((sqrt(e * e - 1.) + (e - 1.) * tan(theta_imp / 2.)) /
                   (sqrt(e * e - 1.) - (e - 1.) * tan(theta_imp / 2.))))
              / k2;
      _t_0 = ((e * sqrt(e * e - 1.) * sin(theta_0)) /
              (1. + e * cos(theta_0)) -
              log((sqrt(e * e - 1.) + (e - 1.) * tan(theta_0 / 2.)) /
                  (sqrt(e * e - 1.) - (e - 1.) * tan(theta_0 / 2.))))
             / k2;
   }
   else
   {
      k2 = sqrt(mu / (8. * pow(r_per, 3.)));

      t_imp = fabs(0.5 * (tan(theta_imp / 2.)
                          + (1. / 3.) * pow(tan(theta_imp / 2.), 3.) / k2));
      _t_0 = fabs(0.5 * (tan(theta_0 / 2.)
                         + (1. / 3.) * pow(tan(theta_0 / 2.), 3.)) / k2);
   }

   vect3dT       r_0v;
   r_0v[0] = -_r_0* cos((pi / 2.) - theta_imp + theta_0);
   r_0v[1] = _r_0 * sin((pi / 2.) - theta_imp + theta_0);
   r_0v[2] = 0.;

   const fType alpha = theta_0 - theta_imp - beta_0;

   vect3dT       v_0v;
   v_0v[0] = -v_0* cos(alpha);
   v_0v[1] = -v_0* sin(alpha);
   v_0v[2] = 0.;

   _r_impv = (_m_tar / m_tot) * r_0v;
   _v_impv = (_m_tar / m_tot) * v_0v;

   _r_tarv = -(_m_imp / m_tot) * r_0v;
   _v_tarv = -(_m_imp / m_tot) * v_0v;
}
}
#endif
