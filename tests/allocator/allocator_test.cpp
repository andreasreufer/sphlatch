#include <iostream>
#include "chunk_allocator.h"
#include "bhtree_nodes.h"

int main()
{
  using namespace sphlatch;

  sleep(1);
  //ChunkAllocator<quadrupoleCellNode, 32768> myAllocator;
  ChunkAllocator<quadrupoleCellNode, 16> myAllocator;
  //ChunkAllocator<int, 8> myAllocator;

  int* intPtr = NULL;

  for (size_t i = 0; i < 32; i++)
  //for (size_t i = 0; i < 100000; i++)
  {
    std::cout << i << "\n";
    myAllocator.pop();
    //intPtr = myAllocator.pop();
    //std::cout << intPtr << "\n";
  }
  sleep(1);
  
  for (size_t i = 0; i < 32; i++)
  //for (size_t i = 0; i < 100000; i++)
  {
    std::cout << i << "\n";
    myAllocator.push();
  }

  sleep(1);
  return 0;
};
