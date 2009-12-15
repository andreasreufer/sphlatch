#include <iostream>
#include <fstream>

#include <boost/lexical_cast.hpp>

//#define SPHLATCH_SINGLEPREC

//#include <omp.h>
//#define SPHLATCH_OPENMP
#define SPHLATCH_HDF5
#define SPHLATCH_NONEIGH

#include "typedefs.h"
typedef sphlatch::fType     fType;
typedef sphlatch::vect3dT   vect3dT;
typedef sphlatch::box3dT    box3dT;

const fType finf = sphlatch::fTypeInf;

#include "impact_solver.cpp"

int main(int argc, char* argv[])
{
#ifdef SPHLATCH_MPI
   MPI::Init(argc, argv);
#endif

   const fType G   = 6.6742e-11;
   const fType m_E = 5.9742e24;

   std::fstream fin;
   fin.open(
      "/Users/areufer/Documents/UniTAPS/collisions/reduced-0-2396511-4.dat",
      std::ios::in);

   std::string token;

   vect3dT r_tar, r_imp, v_tar, v_imp;
   fType   t_0, b_scal, gamma;
   fType   R_tar, R_imp, m_tar, m_imp;

   const fType r_0rel = 10.;

   while (fin)
   {
      fin >> token;

      if (!fin)
         break;

      // header?
      if (token == "#")
      {
         for (size_t j = 1; j <= 45; j++)
            fin >> token;
         continue;
      }

      const size_t ts = boost::lexical_cast<size_t>(token);

      fin >> token;
      const fType st = boost::lexical_cast<fType>(token);

      fin >> token;
      const size_t mergid1 = boost::lexical_cast<size_t>(token);

      fin >> token;
      const size_t mergid2 = boost::lexical_cast<size_t>(token);

      fin >> token;
      const fType mergM = boost::lexical_cast<fType>(token);

      fin >> token;
      const fType mergR = boost::lexical_cast<fType>(token);

      fin >> token;
      const fType mergrho = boost::lexical_cast<fType>(token);

      fin >> token;
      const fType b1m = boost::lexical_cast<fType>(token);

      fin >> token;
      const fType b1R = boost::lexical_cast<fType>(token);

      fin >> token;
      const fType b1rho = boost::lexical_cast<fType>(token);

      fin >> token;
      const size_t b2id1 = boost::lexical_cast<size_t>(token);

      fin >> token;
      const size_t b2id2 = boost::lexical_cast<size_t>(token);

      fin >> token;
      const fType b2m = boost::lexical_cast<fType>(token);

      fin >> token;
      const fType b2R = boost::lexical_cast<fType>(token);

      fin >> token;
      const fType b2rho = boost::lexical_cast<fType>(token);

      fin >> token;
      const fType vrelx = boost::lexical_cast<fType>(token);

      fin >> token;
      const fType vrely = boost::lexical_cast<fType>(token);

      fin >> token;
      const fType vrelz = boost::lexical_cast<fType>(token);

      fin >> token;
      const fType colltime = boost::lexical_cast<fType>(token);

      fin >> token;
      const fType b = boost::lexical_cast<fType>(token);

      fin >> token;
      const fType vimp = boost::lexical_cast<fType>(token);

      const fType v_inf = sqrt(vrelx * vrelx + vrely * vrely + vrelz * vrelz);

      if (b1m > b2m)
      {
         m_tar = b1m;
         m_imp = b2m;
         R_tar = b1R;
         R_imp = b2R;
      }
      else
      {
         m_tar = b2m;
         m_imp = b1m;
         R_tar = b2R;
         R_imp = b1R;
      }

      const fType L_tot = b * m_imp * v_inf;
      const fType r_0   = r_0rel * R_tar;
      const fType gamma = m_imp / m_tar;

      if ((gamma > 0.1) && (m_tar > 0.01 * m_E))
      {
         sphlatch::getImpactSetup(m_tar, m_imp, R_tar, R_imp,
                                  v_inf, L_tot, r_0, G,
                                  r_tar, r_imp, v_tar, v_imp, t_0, b_scal);

         std::cout << m_tar << "\t"
                   << m_imp << "\t"
                   << R_tar << "\t"
                   << R_imp << "\t"
                   << v_inf << "\t"
                   << L_tot << "\t"
                   << r_0 << "\t"
                   << t_0 << "\t"
                   << b   << "\t"
                   << b_scal << "\n";
      }
   }

   fin.close();

#ifdef SPHLATCH_MPI
   MPI::Finalize();
#endif
   return(0);
}
