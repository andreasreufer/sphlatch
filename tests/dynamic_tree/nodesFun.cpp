#include <iostream>
#include <vector>

#include <omp.h>

#include "bhtree.cpp"
typedef sphlatch::BHTree treeType;

typedef sphlatch::particleNode partType;
typedef sphlatch::genericNode nodeType;

typedef sphlatch::quadrupoleCellNode cellType;
typedef sphlatch::quadrupoleCellNode* cellPtr;

typedef sphlatch::costzoneCellNode czllType;
typedef sphlatch::costzoneCellNode* czllPtr;

int main()
{

  typedef nodeType* nodePtrT;
  nodePtrT curPtr;

  std::cout << "instantate new CZ cell ... \n";
  curPtr = new czllType;
  std::cout << "instantate new CZ cell done\n";

  curPtr->isCZ = true;
  curPtr->ident = 42;
  
  std::cout << curPtr->ident << "\n";
  curPtr->clear();
  std::cout << curPtr->ident << "\n";

  static_cast<czllPtr>(curPtr)->cen[0] = 4.2;
  std::cout << static_cast<czllPtr>(curPtr)->cen[0] << " "
            << static_cast<czllPtr>(curPtr)->relCost << "\n";

  static_cast<czllPtr>(curPtr)->clear();

  std::cout << static_cast<czllPtr>(curPtr)->cen[0] << " "
            << static_cast<czllPtr>(curPtr)->relCost << "\n";

  static_cast<czllPtr>(curPtr)->cen[0] = 4.2;
  static_cast<czllPtr>(curPtr)->q23  = 5.3;

  cellType newCell;
  newCell.initFromCZll( * static_cast<czllPtr>(curPtr) );
  
  std::cout << newCell.q23 << "\n";
  
  static_cast<czllPtr>(curPtr)->clear();
  static_cast<czllPtr>(curPtr)->initFromCell(newCell);

  std::cout << static_cast<czllPtr>(curPtr)->cen[0] << " "
            << static_cast<czllPtr>(curPtr)->q23 << "\n";

  std::cout << "    cell size " << sizeof(cellType) << "\n";
  std::cout << " CZ cell size " << sizeof(czllType) << "\n";

  return 0;
}

