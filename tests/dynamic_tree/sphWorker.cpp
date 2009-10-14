#include "sph_fluid_particle.h"
#include "bhtree_particle.h"


namespace sphlatch {

struct SPHfunc
{
   void operator()(treeGhost* _i, const treeGhost* _j)
   {
      _i->pos = _j->pos;
   }
};



template<typename F>
class sphWorker {
public:
   sphWorker() { }
   ~sphWorker() { }

   F sumFunc;

   //void operator()(treeGhost* _i, const treeGhost* _j)
   void operator()()
   {
      sphlatch::treeGhost p1, p2;

      sumFunc(&p1, &p2);
      //F(&p1, &p2);
   }
};
};

int main()
{
   using namespace sphlatch;
   sphWorker<SPHfunc> myWorker;

   myWorker();
}
