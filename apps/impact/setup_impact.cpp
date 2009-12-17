#include <iostream>
#include <fstream>

#include <boost/lexical_cast.hpp>

#include "typedefs.h"
typedef sphlatch::fType     fType;
typedef sphlatch::vect3dT   vect3dT;

#include "impact_solver.cpp"

int main(int argc, char* argv[])
{
#ifdef SPHLATCH_MPI
   MPI::Init(argc, argv);
#endif

   // SI-system: m, s, kg
   const fType G   = 6.6742e-11;
   const fType m_E = 5.9742e24;

   vect3dT r_tarv, r_impv, v_tarv, v_impv;
   fType   t_0, b, b_scal, v_inf, v_imp, e, r_0rel;
   fType   R_tar, R_imp, m_tar, m_imp;

   std::cout << "R_tar [m]  : ";
   std::cin >> R_tar;
   
   std::cout << "R_imp [m]  : ";
   std::cin >> R_imp;

   std::cout << "m_tar [kg] : ";
   std::cin >> m_tar;

   std::cout << "m_imp [kg] : ";
   std::cin >> m_imp;

   std::cout << "v_inf [m/s]: ";
   std::cin >> v_inf;

   std::cout << "b [kgm^2/s]: ";
   std::cin >> b;
   
   std::cout << "r_0 [R_tar]: ";
   std::cin >> r_0rel;

   const fType L_tot = b * m_imp * v_inf;
   const fType r_0   = r_0rel * R_tar;
   const fType gamma = m_imp / m_tar;

   sphlatch::getImpactSetup(m_tar, m_imp, R_tar, R_imp,
                            v_inf, L_tot, r_0, G,
                            r_tarv, r_impv, v_tarv, v_impv, 
                            t_0, b_scal, v_imp, e);

   std::cout << std::scientific;
   std::cout << "m_tar  [kg]     : " << m_tar << "\n"
             << "m_imp  [kg]     : " << m_imp << "\n"
             << "R_tar  [kg]     : " << R_tar << "\n"
             << "R_imp  [kg]     : " << R_imp << "\n"
             << "v_inf  [m/s]    : " << v_inf << "\n"
             << "L_tot  [kgm^2/s]: " << L_tot << "\n"
             << "r_0    [m]      : " << r_0 << "\n"
             << "t_0    [s]      : " << t_0 << "\n"
             << "b      [1]      : " << b << "\n"
             << std::fixed
             << "b_scal [1]      : " << b_scal << "\n";

   std::cout << std::setprecision(3) << std::scientific;
   std::cout << "r_tar0 [m]      : " << r_tarv << "\n"
             << "r_imp0 [m]      : " << r_impv << "\n"
             << "v_tar0 [m/s]    : " << v_tarv << "\n"
             << "v_imp0 [m/s]    : " << v_impv << "\n";

#ifdef SPHLATCH_MPI
   MPI::Finalize();
#endif
   return(0);
}
