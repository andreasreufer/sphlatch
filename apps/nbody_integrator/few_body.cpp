#include <iostream>
#include <sstream>
#include <vector>

#include "typedefs.h"
typedef sphlatch::fType     fType;
typedef sphlatch::cType     cType;
typedef sphlatch::vect3dT   vect3dT;
typedef sphlatch::box3dT    box3dT;

const fType finf = sphlatch::fTypeInf;

///
/// define the particle we are using
///
#include "bhtree_particle.h"
#include "sph_fluid_particle.h"
#include "integrator_predcorr.cpp"
#include "io_particle.h"

class particle :
   public sphlatch::treePart,
   public sphlatch::movingPart,
   public sphlatch::SPHfluidPart,
   public sphlatch::IOPart
{
public:

   sphlatch::PredictorCorrectorO2<vect3dT> posInt;

   void bootstrap()
   {
      posInt.bootstrap(pos, vel, acc);
   }

   void predict(const fType _dt)
   {
      posInt.predict(pos, vel, acc, _dt);
   }

   void correct(const fType _dt)
   {
      posInt.correct(pos, vel, acc, _dt);
   }

   ioVarLT getLoadVars() { ioVarLT vars;
                           return(vars); }
   ioVarLT getSaveVars() { ioVarLT vars;
                           return(vars); }
};

typedef particle                       partT;

#include "particle_set.cpp"
typedef sphlatch::ParticleSet<partT>   partSetT;

// particles are global
partSetT parts;

void derive()
{
   const int   nop = parts.getNop();
   const fType G   = parts.attributes["gravconst"];

   vect3dT cacc;

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
            //const fType h = 0.5*(parts[i].h + parts[j].h );
            const fType mOr3 = parts[i].m / (r * r * r);

            cacc -= mOr3 * rvec;
         }
      }

      parts[j].acc = G * cacc;
   }
}

fType timestep(const fType _stepTime)
{
   const fType time = parts.attributes["time"];

   fType dtSave = (floor((time / _stepTime) + 1.e-6) + 1.) * _stepTime - time;
   fType dtA    = finf;

   const size_t nop = parts.getNop();

   for (size_t i = 0; i < nop; i++)
   {
      const fType ai = sqrt(dot(parts[i].acc, parts[i].acc));
      if (ai > 0.)
      {
         const fType dtAi = sqrt(parts[i].h / ai);
         dtA = dtAi < dtA ? dtAi : dtA;
      }
   }

   dtA *= 0.5;
   dtA = 10.;

   //FIXME: globally minimize all dts
   fType dt = finf;

   dt = dtSave < dt ? dtSave : dt;
   dt = dtA < dt ? dtA : dt;
   return(dt);
}

int main(int argc, char* argv[])
{
   if (not (argc == 3))
   {
      std::cerr << "usage: few_body <outStepTime> <stopTime>\n";
      return(1);
   }

   std::istringstream stepStr(argv[1]);
   fType stepTime;
   stepStr >> stepTime;

   std::istringstream stopStr(argv[2]);
   fType stopTime;
   stopStr >> stopTime;

   parts.resize(2);


   fType& time(parts.attributes["time"]);
   cType& step(parts.step);

   parts.attributes["time"]      = -3200;
   parts.attributes["gravconst"] = 6.6742e-8;

   parts[0].pos = 1.084204e+09, -2.036850e+04, 2.400897e+09;
   parts[0].vel = -3.029064e+05, -4.146875e+00, -3.616368e+05;
   parts[0].m   = 8.584803e+26;

   parts[1].pos = -1.617881e+08, -2.427607e+04, -3.582495e+08;
   parts[1].vel = 4.520210e+04, 6.975465e-01, 5.396420e+04;
   parts[1].m   = 5.752381e+27;

   const size_t nop = parts.getNop();

   // bootstrapping the integrator
   derive();
   for (size_t i = 0; i < nop; i++)
      parts[i].bootstrap();

   // start the loop
   while (time < stopTime)
   {
      derive();
      const fType dt = timestep(stepTime);
      for (size_t i = 0; i < nop; i++)
         parts[i].predict(dt);

      time += dt;

      derive();
      for (size_t i = 0; i < nop; i++)
         parts[i].correct(dt);
      step++;

      std::cerr << "corrected (t = " << time << ")\n";

      const fType dumpTime = (floor((time / stepTime) + 1.e-9)) * stepTime;
      if (fabs(dumpTime - time) < 1.e-9)
      {
         std::cout << time << "\t";
         for (size_t i = 0; i < nop; i++)
         {
            std::cout << std::scientific;
            std::cout << parts[i].pos[0] << "\t"
                      << parts[i].pos[1] << "\t"
                      << parts[i].pos[2] << "\t"
                      << parts[i].vel[0] << "\t"
                      << parts[i].vel[1] << "\t"
                      << parts[i].vel[2] << "\t"
                      << parts[i].acc[0] << "\t"
                      << parts[i].acc[1] << "\t"
                      << parts[i].acc[2] << "\t";
         }
         std::cout << "\n";
      }
   }

   return(0);
}
