#include <iostream>
#include <vector>

#include <omp.h>

#define SPHLATCH_PADD64
//include "bhtree_dynamic.h"
//#include "bhtree_nodes.h"
//#include "bhtree_particle.h"
//typedef sphlatch::BHTree               treeType;
#include "bhtree.cpp"

typedef sphlatch::genericNode          nodeType;
typedef sphlatch::particleNode         pnodType;
typedef sphlatch::genericCellNode      gcllType;
typedef sphlatch::monopoleCellNode     monoType;
typedef sphlatch::quadrupoleCellNode   cellType;
typedef sphlatch::costzoneCellNode     czllType;

typedef sphlatch::treePart         partType;

int main()
{
   std::cout << " generic   node " << sizeof(nodeType) << "\n"
             << " particle  node " << sizeof(pnodType) << "\n"
             << " gen  cell node " << sizeof(gcllType) << "\n"
             << " mono cell node " << sizeof(monoType) << "\n"
             << " quad cell node " << sizeof(cellType) << "\n"
             << " CZ   cell node " << sizeof(czllType) << "\n"
             << " particle       " << sizeof(partType) << "\n";

   return(0);
}
