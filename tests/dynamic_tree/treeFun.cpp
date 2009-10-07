#include <iostream>
#include <vector>

#include <omp.h>

#include "bhtree.cpp"
typedef sphlatch::BHTree treeT;

typedef sphlatch::nodeT       nodeT;
typedef sphlatch::pnodT       pnodT;
typedef sphlatch::gcllT       gcllT;
typedef sphlatch::mcllT       mcllT;
typedef sphlatch::qcllT       qcllT;
typedef sphlatch::czllT       czllT;

typedef sphlatch::treeGhost   partT;

/*
#include "bhtree_part_insertmover.h"
typedef sphlatch::BHTreePartsInsertMover inserterT;

#include "bhtree_cz_builder.h"

#include "bhtree_worker_grav.h"
typedef sphlatch::BHTreeWorkerGrav gravworkerT;

#include "bhtree_treedump.h"
typedef sphlatch::BHTreeDump dumpT;*/

#include "io_manager.h"
typedef sphlatch::IOManager io_type;

#include "particle_manager.h"
typedef sphlatch::ParticleManager part_type;

int main()
{
  part_type& PartManager(part_type::instance());
  io_type& IOManager(io_type::instance());

  /*
  sleep(1);
  treeT& Tree(treeT::instance());
  sleep(1);

  IOManager.loadDump("test.h5part");

  dumpT dumper(&Tree);
  dumper.dotDump("test.dot");

  inserterT inserter(&Tree);
  
  gravworkerT worker(&Tree);
#pragma omp parallel for firstprivate(worker)
  for (int i = 0; i < 8; i++)
  {
    const size_t tid = omp_get_thread_num();
    std::cout << tid << ":" << i << ":" << &worker << "\n";
  }
  sleep(1);*/
  return 0;
}

