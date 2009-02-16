#include <iostream>
#include <vector>

#include <omp.h>

#define SPHLATCH_PADTO64BYTES
#include "bhtree_dynamic.h"
typedef sphlatch::BHTree               treeType;

typedef sphlatch::genericNode          nodeType;
typedef sphlatch::particleNode         partType;
typedef sphlatch::genericCellNode      gcllType;
typedef sphlatch::monopoleCellNode     monoType;
typedef sphlatch::quadrupoleCellNode   cellType;
typedef sphlatch::costzoneCellNode     czllType;

int main()
{
   std::cout << " generic   node " << sizeof(nodeType) << "\n"
             << " particle  node " << sizeof(partType) << "\n"
             << " gen  cell node " << sizeof(gcllType) << "\n"
             << " mono cell node " << sizeof(monoType) << "\n"
             << " quad cell node " << sizeof(cellType) << "\n"
             << " CZ   cell node " << sizeof(czllType) << "\n";

   return(0);
}
