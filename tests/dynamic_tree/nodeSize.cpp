#include <iostream>
#include <vector>

#include <omp.h>

#define SPHLATCH_PADD64
//#include "bhtree.cpp"

#include "typedefs.h"

//#include "bhtree_nodes.cpp"
#include "bhtree.cpp"

//typedef sphlatch::BHTree               treeT;

typedef sphlatch::nodeT nodeT;
typedef sphlatch::pnodT pnodT;
typedef sphlatch::gcllT gcllT;
typedef sphlatch::mcllT mcllT;
typedef sphlatch::qcllT qcllT;
typedef sphlatch::czllT czllT;

typedef sphlatch::treeGhost        partT;

int main()
{
   std::cout << " generic   node " << sizeof(nodeT) << "\n"
             << " particle  node " << sizeof(pnodT) << "\n"
             << " gen  cell node " << sizeof(gcllT) << "\n"
             << " mono cell node " << sizeof(mcllT) << "\n"
             << " quad cell node " << sizeof(qcllT) << "\n"
             << " CZ   cell node " << sizeof(czllT) << "\n"
             << " particle       " << sizeof(partT) << "\n";

   return(0);
}
