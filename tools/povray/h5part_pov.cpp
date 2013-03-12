#include <iostream>
#include <sstream>
#include <vector>

//#define SPHLATCH_SINGLEPREC

#include <omp.h>
#define SPHLATCH_OPENMP
#define SPHLATCH_HDF5
#define SPHLATCH_NONEIGH

#include "typedefs.h"
typedef sphlatch::fType     fType;
typedef sphlatch::cType     cType;
typedef sphlatch::vect3dT   vect3dT;
typedef sphlatch::box3dT    box3dT;
typedef sphlatch::partsIndexListT   plistT;


const fType finf = sphlatch::fTypeInf;

#include "bhtree.cpp"
typedef sphlatch::BHTree    treeT;

///
/// define the particle we are using
///
#include "bhtree_particle.h"
#include "sph_fluid_particle.h"
#include "io_particle.h"

class particle :
   public sphlatch::treePart,
   public sphlatch::SPHfluidPart,
   public sphlatch::IOPart,
   public sphlatch::ANEOSPart
{
public:
   ioVarLT getLoadVars()
   {
      ioVarLT vars;
      
      vars.push_back(storeVar(id,  "id"));
      vars.push_back(storeVar(pos, "pos"));
      vars.push_back(storeVar(mat, "mat"));
      vars.push_back(storeVar(h, "h"));

      return(vars);
   }

   ioVarLT getSaveVars()
   {
      ioVarLT vars;
      return(vars);
   }
};

typedef particle                       partT;

#include "particle_set.cpp"
typedef sphlatch::ParticleSet<partT>   partSetT;

// particles are global
partSetT parts;

//#include "misc_physics.cpp"
//using sphlatch::addToCOM;

int main(int argc, char* argv[])
{
#ifdef SPHLATCH_MPI
   MPI::Init(argc, argv);
#endif

   if (not (argc == 5))
   {
      std::cerr <<
      "usage: cut_sphere <input dump> <output dump> <rmax> <vmax>\n";
      return(1);
   }


   std::string iname = argv[1];
   std::string oname = argv[2];

   std::istringstream rmaxStr(argv[3]);
   fType rmax;
   rmaxStr >> rmax;

   std::istringstream vmaxStr(argv[4]);
   fType vmax;
   vmaxStr >> vmax;

   // load the particles
   parts.loadHDF5(iname);
   const size_t nop = parts.getNop();

   std::cerr << "loaded " << nop << " particles\n";

   std::fstream fout;
   fout.open(oname.c_str(), std::ios::out);

   if (!fout)
     return(1);

   const fType scl = 1.e7;
   //const fType r = 9520359.;
   const fType r = 0.3*9520359./scl;

   fout << "#version 3.6;\n";
   fout << "global_settings {  assumed_gamma 1.0 }\n";
   fout << "#default{ finish{ ambient 0.1 diffuse 0.9}}\n";
   fout << "#include \"colors.inc\"\n";
   fout << "#declare tsilc = rgb <0.54509804,0.,0.>;\n";
   fout << "#declare twatr = rgb <0.,0.54509804,0.54509804>;\n";
   fout << "#declare tiron = rgb <0.,0.,0.80392157>;\n";
   fout << "#declare isilc = rgb <1.,0.54901961,0.>;\n";
   fout << "#declare iwatr = rgb <0.,1.,1.>;\n";
   fout << "#declare iiron = rgb <0.11764706,0.56470588,1.>;\n";
   fout << "camera {location < 0.0, 0.0,80.0 >\n";
   fout << "        up       < 0.0, 4.5, 0.0 > \n";
   fout << "        right    < 8.0, 0.0, 0.0 > \n";
   fout << "        look_at  < 0.0, 0.0, 0.0 > \n";
   fout << "        rotate   <   0,   0,  90> } \n";
   fout << "light_source{<1500,2500,  500> color White}\n";
   //fout << "light_source{<1.5e11, 2.5e11, -2.5e11> color White}\n";

   std::vector<std::string> colors;

   colors.resize(32);

   colors[1] = "tsilc";
   colors[2] = "twatr";
   colors[5] = "tiron";
   
   colors[1 + 16] = "isilc";
   colors[2 + 16] = "iwatr";
   colors[5 + 16] = "iiron";

   for (size_t i = 0; i < nop; i++)
   {
     size_t coli = parts[i].mat;
     if (parts[i].id > 2.e6)
       coli += 16;

     if (parts[i].pos[2] < 0.)
     {
     fout << "sphere{<" 
          << parts[i].pos[0] / scl << ","
          << parts[i].pos[1] / scl << ","
          << parts[i].pos[2] / scl << ">, "
          << r << "\n"
          << "pigment {color " << colors[coli] << "}\n} \n";
     }
   }

   fout.close();

   return(0);
}
