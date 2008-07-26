#ifndef SPHLATCH_EOS_TILLOTSON
#define SPHLATCH_EOS_TILLOTSON

#include "typedefs.h"
#include "eos_generic.h"

namespace sphlatch {

class Tillotson : public EOS<Tillotson> {

public:
void pressure()
{
  std::cout << "pressure Tillotson!\n";
};

void init()
{
  Logger << "init Tillotson EOS";
};

};
}
#endif
