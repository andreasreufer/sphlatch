#ifndef SPHLATCH_EOS_IDEALGAS
#define SPHLATCH_EOS_IDEALGAS

#include "typedefs.h"
#include "eos_generic.h"

namespace sphlatch {

class IdealGas : public EOS<IdealGas> {

public:
void pressure()
{
  std::cout << "pressure!\n";
};

void init()
{
};

};
}
#endif
