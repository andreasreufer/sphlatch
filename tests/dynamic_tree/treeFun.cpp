#include <iostream>
#include <vector>

#include <omp.h>

#include "bhtree_dynamic.h"
#include "bhtree_worker.h"
#include "bhtree_worker_grav.h"
#include "bhtree_part_inserter.h"
//#include "bhtree_cz_builder.h"
#include "bhtree_worker_treedump.h"

typedef sphlatch::particleNode partType;
typedef sphlatch::genericNode nodeType;
typedef sphlatch::quadrupoleCellNode cellType;
typedef sphlatch::quadrupoleCellNode* cellPtr;
typedef sphlatch::costzoneCellNode czllType;

typedef sphlatch::BHTree treeType;
typedef sphlatch::BHTreeWorker workerType;
typedef sphlatch::BHTreeDump dumpType;
typedef sphlatch::BHTreeWorkerGrav gravworkerType;
//typedef sphlatch::BHTreePartsInserter inserterType;

int main()
{
  treeType& Tree(treeType::instance());

  //inserterType inserter;

  dumpType dumper(&Tree);
  dumper.dotDump("test.dot");

  //workerType worker(&Tree);
  gravworkerType worker(&Tree);
#pragma omp parallel for firstprivate(worker)
  for (int i = 0; i < 8; i++)
  {
    //const size_t tid = omp_get_thread_num();
    //std::cout << tid << ":" << i << ":" << &worker << "\n";
    
  }
  return 0;
}

