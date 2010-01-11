#include <iostream>
#include <vector>

#define SPHLATCH_HDF5
#define SPHLATCH_NONEIGH

#include "typedefs.h"
typedef sphlatch::fType     fType;
typedef sphlatch::vect3dT   vect3dT;

///
/// define the particle we are using
///
#include "bhtree_particle.h"
#include "sph_fluid_particle.h"
#include "io_particle.h"

class particle :
   public sphlatch::treePart,
   public sphlatch::movingPart,
   public sphlatch::SPHfluidPart,
   public sphlatch::IOPart
{
public:

   ioVarLT getLoadVars()
   {
      ioVarLT vars;

      vars.push_back(storeVar(pos, "pos"));
      vars.push_back(storeVar(m, "m"));
      vars.push_back(storeVar(h, "h"));
      vars.push_back(storeVar(id, "id"));
      return(vars);
   }

   ioVarLT getSaveVars()
   {
      ioVarLT vars;

      vars.push_back(storeVar(pos, "pos"));
      vars.push_back(storeVar(acc, "acc"));
      vars.push_back(storeVar(m, "m"));
      vars.push_back(storeVar(h, "h"));

      vars.push_back(storeVar(id, "id"));

      return(vars);
   }
};

typedef particle                       partT;

#include "particle_set.cpp"
typedef sphlatch::ParticleSet<partT>   partSetT;

// particles are global
partSetT parts;

fType splineOSmoR3(const fType _r, const fType _h)
{
   const fType u = _r / _h;

   if (u >= 2.)
   {
      const fType r3 = _r * _r * _r;
      return(1. / r3);
   }
   else
   if (u > 1.)
   {
      const fType r3 = _r * _r * _r;
      return((1. / r3) * (
                -(1. / 15.)
                + (8. / 3.) * u * u * u
                - 3. * u * u * u * u
                + (6. / 5.) * u * u * u * u * u
                - (1. / 6.) * u * u * u * u * u * u
                ));
   }
   else
   {
      const fType h3 = _h * _h * _h;
      return((1. / h3) * (
                +(4. / 3.)
                - (6. / 5.) * u * u
                + (1. / 2.) * u * u * u
                ));
   }
}

void derive()
{
   const int   nop = parts.getNop();
   const fType G   = parts.attributes["gravconst"];

   vect3dT cacc;

#pragma omp parallel for private(cacc)
   for (int j = 0; j < nop; j++)
   {
      cacc = 0., 0., 0.;

      for (int i = 0; i < nop; i++)
      {
         if (i != j)
         {
            const vect3dT rvec = parts[j].pos - parts[i].pos;
            const fType   r    = sqrt(rvec[0] * rvec[0] +
                                      rvec[1] * rvec[1] +
                                      rvec[2] * rvec[2]);
            const fType h    = parts[i].h;
            const fType mOr3 = parts[i].m* splineOSmoR3(r, h);

            cacc -= mOr3 * rvec;
         }
      }

      parts[j].acc = G * cacc;
   }
}

int main(int argc, char* argv[])
{
   if (not argc == 3)
   {
      std::cerr << "usage: nbody_bf <inputdump> <outputdump>\n";
      return(1);
   }

   std::string inFilename  = argv[1];
   std::string outFilename = argv[2];


   // load the particles
   parts.loadHDF5(inFilename);

   // derive
   derive();

   // store the particles
   parts.doublePrecOut();
   parts.saveHDF5(outFilename);

   return(0);
}
