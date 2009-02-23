#include <iostream>
#include <vector>

#include <omp.h>

#define SPHLATCH_PADTO64BYTES
#include "bhtree_dynamic.h"
typedef sphlatch::BHTree               treeType;

typedef sphlatch::genericNode          nodeType;
typedef sphlatch::particleNode         pnodType;
typedef sphlatch::genericCellNode      gcllType;
typedef sphlatch::monopoleCellNode     monoType;
typedef sphlatch::quadrupoleCellNode   cellType;
typedef sphlatch::costzoneCellNode     czllType;

typedef sphlatch::particleGeneric      partType;

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
