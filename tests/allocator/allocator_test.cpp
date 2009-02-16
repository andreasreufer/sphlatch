#include <iostream>

#include "allocator_chunk.h"
#include "allocator_simple.h"
#include "allocator_fixedbucket.h"
#include "bhtree_node_cells.h"

int main()
{
  using namespace sphlatch;

  sleep(1);
  //ChunkAllocator<quadrupoleCellNode, 32768> myAllocator;
  //ChunkAllocator<quadrupoleCellNode, 16> myAllocator;
  ChunkAllocator<int, 65535> myAllocator;
  //SimpleAllocator<int, 65535> myAllocator;

  int* intPtr = NULL;

  for (size_t i = 0; i < 32; i++)
  //for (size_t i = 0; i < 100000; i++)
  {
    intPtr = myAllocator.pop();
    *intPtr = i;
    std::cout << i << ": " << *intPtr << "\n";
    //intPtr = myAllocator.pop();
    //std::cout << intPtr << "\n";
  }
  sleep(1);
  
  for (size_t i = 0; i < 32; i++)
  //for (size_t i = 0; i < 100000; i++)
  {
    std::cout << i << "\n";
    //myAllocator.push(intPtr);
  }

  sleep(1);
  return 0;
};
