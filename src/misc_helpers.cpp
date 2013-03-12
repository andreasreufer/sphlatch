#ifndef SPHLATCH_MISC_HELPERS
#define SPHLATCH_MISC_HELPERS

#include <sys/stat.h>

namespace sphlatch
{

  bool fileExists(std::string _fname)
  {
    struct stat statBuff;

    if ( stat(_fname.c_str(), &statBuff) == -1)
      return false;
    else
      return true;
  }
}

#endif
