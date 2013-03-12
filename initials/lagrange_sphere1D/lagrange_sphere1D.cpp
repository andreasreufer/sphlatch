// uncomment for single-precision calculation
//#define SPHLATCH_SINGLEPREC

// enable parallel version
//#define SPHLATCH_PARALLEL

#define SPHLATCH_HDF5
#define SPHLATCH_LOGGER


#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

#include <boost/lexical_cast.hpp>

#include "typedefs.h"
typedef sphlatch::fType                    fType;
typedef sphlatch::iType                    iType;

typedef sphlatch::fvectT                   fvectT;
typedef sphlatch::idvectType               idvectT;

#include "constants.h"

#include "hdf5_io.cpp"
typedef sphlatch::HDF5File                 hdf5fT;

#include "lagrange_sphere1D_solver.cpp"
typedef sphlatch::LagrangeSphere1DSolver   lg1D_solverT;

using namespace sphlatch::constants;
using namespace sphlatch::constants::unitsCGS;

int main(int argc, char* argv[])
{
   if ( argc != 5 )
   {
     std::cerr << "lagrange_sphere1D <infile> <outfile> <stoptime> <fric>\n";
     return(EXIT_FAILURE);
   }


   std::istringstream stopStr(argv[3]);
   fType stopTime;
   stopStr >> stopTime;
   
   std::istringstream fricStr(argv[4]);
   fType fricTime;
   fricStr >> fricTime;

   hdf5fT InFile(argv[1]);

   size_t noCells, idummy;
   InFile.getDims("m", noCells, idummy);
   const size_t noEdges = noCells + 1;
   
   ///
   /// instantate solver
   ///
   lg1D_solverT Solver(noCells);

   // edge values
   fvectT& r(Solver.r);
   fvectT& v(Solver.v);
   
   // centered values
   fvectT& m(Solver.m);
   fvectT& rho(Solver.rho);
   fvectT& u(Solver.u);
   fvectT& p(Solver.p);
   fvectT& dudt(Solver.dudt);
   fvectT& T(Solver.T);
   fvectT& S(Solver.S);
   fvectT& cs(Solver.cs);

   idvectT& mat(Solver.mat);
   idvectT& phase(Solver.phase);

   InFile.loadPrimitive("m", m);
   InFile.loadPrimitive("rho", rho);
   InFile.loadPrimitive("S", S);
   InFile.loadPrimitive("mat", mat);
   InFile.loadPrimitive("p", p);
   Solver.gravConst = InFile.loadAttribute("gravconst");

   ///
   /// set the shell edges
   ///
   r(0) = 0.;
   v(0) = 0.;
   for (size_t i = 1; i < noEdges; i++)
   {
      r(i) = pow(pow(r(i-1), 3.) +
                 (3. / (4 * M_PI)) * (m(i-1) / rho(i-1)), 1. / 3.);
      v(i) = 0.;
   }
   std::cout << r(1) << "\n";

   ///
   /// integrate for a certain physical time
   ///
   std::cerr << " start 1D Lagrange solver\n";
   Solver.friction = 1./fricTime;
   Solver.integrateTo(stopTime);
   std::cerr << " ... finished\n";
   
   ///
   /// the last cell contains vacuum, which will assign
   /// far outside zero mass. circumvent this, by setting
   /// the same density and temperature like on the second
   /// to last cell
   ///
   rho(noCells - 1) = rho(noCells - 2);
   u(noCells - 1)   = u(noCells - 2);

   ///
   /// calculate the position of the cell centers
   ///
   fvectT rc, vc;
   rc.resize(noCells);
   vc.resize(noCells);

   for (size_t i = 0; i < noCells; i++)
   {
      rc(i) = 0.5 * (r(i) + r(i + 1));
      vc(i) = 0.5 * (v(i) + v(i + 1));
   }

   hdf5fT OutFile(argv[2]);
   OutFile.savePrimitive("r", rc);
   OutFile.savePrimitive("v", vc);
   OutFile.savePrimitive("m", m);
   OutFile.savePrimitive("u", u);
   OutFile.savePrimitive("p", p);
   OutFile.savePrimitive("S", S);
   OutFile.savePrimitive("mat", mat);
   OutFile.savePrimitive("rho", rho);
   OutFile.savePrimitive("dudt", dudt);
   OutFile.savePrimitive("T", T);
   OutFile.savePrimitive("phase", phase);
   
   OutFile.saveAttribute("gravconst", Solver.gravConst);
   OutFile.saveAttribute("frictime", 1./Solver.friction);

   return(EXIT_SUCCESS);
}
