#include <iostream>
#include <vector>

#include <omp.h>

#define SPHLATCH_PADD64

#include "bhtree.cpp"
#include "typedefs.h"

typedef sphlatch::nodeT       nodeT;
typedef sphlatch::pnodT       pnodT;
typedef sphlatch::gcllT       gcllT;
typedef sphlatch::mcllT       mcllT;
typedef sphlatch::qcllT       qcllT;
typedef sphlatch::czllT       czllT;

typedef sphlatch::treeGhost   partT;

typedef sphlatch::BHTree      treeT;

int main()
{
   std::cout << " generic   node " << sizeof(nodeT) << "   "
             << (sizeof(nodeT) % 64 ) << "\n"
             << " particle  node " << sizeof(pnodT) << "   "
             << (sizeof(pnodT) % 64 ) << "\n"
             << " gen  cell node " << sizeof(gcllT) << "   " 
             << (sizeof(gcllT) % 64 ) << "\n"
             << " mono cell node " << sizeof(mcllT) << "   " 
             << (sizeof(mcllT) % 64 ) << "\n"
             << " quad cell node " << sizeof(qcllT) << "   " 
             << (sizeof(qcllT) % 64 ) << "\n"
             << " CZ   cell node " << sizeof(czllT) << "   " 
             << (sizeof(czllT) % 64 ) << "\n"
             << " particle       " << sizeof(partT) << "   " 
             << (sizeof(partT) % 64 ) << "\n";

   return(0);
}
