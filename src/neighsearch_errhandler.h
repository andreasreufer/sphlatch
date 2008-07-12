#ifndef SPHLATCH_NEIGHSEARCH_ERRHANDLER_H
#define SPHLATCH_NEIGHSEARCH_ERRHANDLER_H

/*
 *  neighsearch_errhandler.h
 *
 *
 *  Created by Andreas Reufer on 10.07.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"
#include "err_handler.h"
#include <boost/lexical_cast.hpp>

namespace sphlatch
{
class TooManyNeighs : public GenericError
{
public:
typedef genericNode* nodePtrT;
typedef genericCellNode* cellPtrT;

TooManyNeighs(size_t _part,
              size_t _noNeighs,
              partsIndexVectType _neighList,
              valvectRefType _neighDist);

~TooManyNeighs();
private:
  idvectType   neighID, partID;
  valvectType  neighDist;
  matrixType neighPos, partPos;

  std::string dumpFile;
};

TooManyNeighs::TooManyNeighs(size_t _part,
                             size_t _noNeighs,
                             partsIndexVectType _neighList,
                             valvectRefType _neighDist)
{
  using namespace sphlatch::vectindices;
  idvectRefType id(PartManager.id);
  matrixRefType pos(PartManager.pos);

  neighID.resize(_noNeighs);
  partID.resize(1);

  neighDist.resize(_noNeighs);
  neighPos.resize(_noNeighs, 3);
  partPos.resize(1,3);

  dumpFile.clear();
  dumpFile.append("neighbours_part");
  dumpFile.append( boost::lexical_cast<std::string>( id(_part) ) );
  dumpFile.append(".hdf5");

  partID(0) = id(_part);
  partPos(0, X) = pos(_part, X);
  partPos(0, Y) = pos(_part, Y);
  partPos(0, Z) = pos(_part, Z);

  for (size_t i = 0; i < _noNeighs; i++)
  {
    const size_t curNeigh = _neighList[i+1];
    neighID(i) = id( curNeigh );
    neighDist(i) = _neighDist(i+1);

    neighPos(i, X) = pos( curNeigh, X);
    neighPos(i, Y) = pos( curNeigh, Y);
    neighPos(i, Z) = pos( curNeigh, Z);
  }

  IOManager.savePrimitive( partID,    "part_id",   dumpFile);
  IOManager.savePrimitive( partPos,   "part_pos",  dumpFile);
  IOManager.savePrimitive( neighID,   "neighs_id",   dumpFile);
  IOManager.savePrimitive( neighPos,  "neighs_pos",  dumpFile);
  IOManager.savePrimitive( neighDist, "neighs_dist", dumpFile);

  Logger.stream << "error: particle " << id(_part)
                << " has too many neighbours ("
                << _noNeighs << "), stored them to "
                << dumpFile;
  Logger.flushStream();
  Logger.destroy();
};

TooManyNeighs::~TooManyNeighs()
{
};
};

#endif
