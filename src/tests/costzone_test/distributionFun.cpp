#include <iostream>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include "typedefs.h"

int main()
{
  const size_t noCells = 32768;

  std::string readToken;

  sphlatch::countsVectType partsPerCell;

  partsPerCell.resize(noCells);

  std::fstream fin;
  fin.open( "cells.txt", std::ios::in);
  
  for (size_t i = 0; i < noCells; i++)
  {
    fin >> readToken;
    fin >> readToken;
    partsPerCell[i] = boost::lexical_cast<size_t>(readToken);
    std::cout << partsPerCell[i] << "\n";
    
    fin >> readToken;
  }

  return 0;
};
