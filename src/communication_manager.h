#ifndef SPHLATCH_COMMUNICATION_MANAGER_H
#define SPHLATCH_COMMUNICATION_MANAGER_H

#include <vector>
#include <list>
#include <mpi.h>

#include <boost/assign/std/vector.hpp>

#include "particle_manager.h"

namespace sphlatch
{
class CommunicationManager
{
public:
typedef CommunicationManager self_type;
typedef CommunicationManager& self_reference;
typedef CommunicationManager* self_pointer;

///
/// send queue element: vector iterators for partIndex vector,
/// number of particles to send and rank where to send to
///
struct sendQueueElemType {
  partsIndexVectType::const_iterator itr;
  partsIndexVectType::const_iterator end;
  size_t noParts;
  size_t sendRank;
};

///
/// we will use a list
///
typedef std::list< sendQueueElemType > sendQueueType;
typedef sendQueueType& sendQueueRefType;

typedef sphlatch::ParticleManager PartManagerType;

static self_reference instance(void);

PartManagerType& PartManager;

///
/// this is the central exchange part
///
void exchange(domainPartsIndexRefType _partsIndices,
              size_t _noGhosts,
              quantsTypeRefType _quantities);

void sendGhostsPrepare(domainPartsIndexRefType _ghostsIndices);

void sendGhosts(idvectRefType _idVect);
void sendGhosts(valvectRefType _valVect);
void sendGhosts(matrixRefType _matrix);

void sendGhosts(quantsTypeRefType _quantities);

private:
std::vector<MPI::Request> recvReqs;
template<class T> void queuedExch(T& _src, T& _trgt,
                                  const sendQueueRefType _queue,
                                  countsVectType _offsets,
                                  countsVectType _noRecvParts);

void sendChunk(const idvectRefType _src, sendQueueElemType _qelem);
idvectType idSendBuff;
void sendChunk(const valvectRefType _src, sendQueueElemType _qelem);
valvectType scalSendBuff;
void sendChunk(const matrixRefType _src, sendQueueElemType _qelem);
matrixType vectSendBuff;

void recvChunk(idvectRefType _target, countsRefType _offset,
               MPI::Request& _recvRq,
               const size_t& _noParts, const size_t& _recvFrom);
void recvChunk(valvectRefType _target, countsRefType _offset,
               MPI::Request& _recvRq,
               const size_t& _noParts, const size_t& _recvFrom);
void recvChunk(matrixRefType _target, countsRefType _offset,
               MPI::Request& _recvRq,
               const size_t& _noParts, const size_t& _recvFrom);

void prepareQueuesOffsets(domainPartsIndexRefType _indices,
                          size_t _totalOffset,
                          sendQueueRefType _queue,
                          countsVectRefType _offsets,
                          countsVectRefType _noRecvParts);

sendQueueType ghostQueue;
countsVectType localOffsets, noRecvParts, ghostOffsets, noRecvGhosts;
size_t identTypeSize, valueTypeSize, newNoParts, newNoGhosts;

public:
///
/// various simple helper functions
///
void sendMatrix(matrixRefType _matrix, size_t _recvDomain);
void recvMatrix(matrixRefType _matrix, size_t _sendDomain);

void sendBitset(bitsetRefType _bitset, size_t _recvDomain);
void recvBitset(bitsetRefType _bitset, size_t _sendDomain);

template<class T> void sendVector(std::vector<T>& _vector, size_t _recvDomain);
template<class T> void recvVector(std::vector<T>& _vector, size_t _sendDomain);

void sumUpCounts(countsVectRefType _indexVect);

void max(valueRefType _val);
void min(valueRefType _val);
void sum(valueRefType _val);

size_t getMyDomain();
size_t getNoDomains();

void barrier();

size_t domainToMPIrank(size_t _domain);
size_t MPIrankToDomain(size_t _rank);

protected:
CommunicationManager(void);
~CommunicationManager(void);

private:
static self_pointer _instance;
static size_t commBuffSize, noBuffParts;
size_t myDomain, myRank, noDomains;

std::vector<size_t> domainToRank;
std::vector<size_t> rankToDomain;
};

CommunicationManager::self_pointer CommunicationManager::_instance = NULL;

CommunicationManager::self_reference CommunicationManager::instance(void)
{
  if (_instance != NULL)
    return *_instance;
  else
    {
      _instance = new CommunicationManager;
      return *_instance;
    }
}

CommunicationManager::CommunicationManager(void) :
  PartManager(PartManagerType::instance())
{
  noDomains = MPI::COMM_WORLD.Get_size();
  domainToRank.resize(noDomains);
  rankToDomain.resize(noDomains);

  ///
  /// insert mapping here
  ///
  for (size_t i = 0; i < noDomains; i++)
    {
      domainToRank[i] = i;
      rankToDomain[i] = i;
    }

  myRank = MPI::COMM_WORLD.Get_rank();
  myDomain = rankToDomain[ myRank ];

  localOffsets.resize(noDomains);
  ghostOffsets.resize(noDomains);

  noRecvParts.resize(noDomains);
  noRecvGhosts.resize(noDomains);

  recvReqs.resize(noDomains);

  identTypeSize = sizeof(identType);
  valueTypeSize = sizeof(valueType);
}

CommunicationManager::~CommunicationManager(void)
{
}

///
/// buffer size in bytes, MPI transfers will have on average this size
///
//size_t CommunicationManager::commBuffSize =  131072;  //  128 kByte
size_t CommunicationManager::commBuffSize = 1048576;    // 1024 kByte
//size_t CommunicationManager::commBuffSize = 8388608;  // 8096 kByte

///
/// buffer size in particles for the exchange() functions
/// 64k gives about 1.5MB for a 3D double precision vector
///
/// you may want to lower this value, when you plan to transmit
/// higher dimension vectors like tensors
///
//size_t CommunicationManager::noBuffParts =   200; // for testing purposes
size_t CommunicationManager::noBuffParts = 65536;

void CommunicationManager::exchange(domainPartsIndexRefType _partsIndices,
                                    size_t _noGhosts,
                                    quantsTypeRefType _quantities)
{
  ///
  /// prepare particle queue and number of particles
  ///
  sendQueueType partsQueue;

  prepareQueuesOffsets(_partsIndices, 0, partsQueue,
                       localOffsets, noRecvParts);

  newNoParts = 0;
  for (size_t i = 0; i < noDomains; i++)
    {
      newNoParts += noRecvParts[i];
    }

  rangeType newLocRange(0, newNoParts);

  ///
  /// determine information for particles staying on node
  ///
  const size_t stayOffset = localOffsets[myRank];
  noRecvParts[myRank] = 0;

  partsIndexVectType::const_iterator stayItr;
  partsIndexVectType::const_iterator stayEnd = _partsIndices[myDomain].end();

  ///
  /// set new number of particles (noGhosts needs to be known)
  ///
  PartManager.setNoParts(newNoParts, _noGhosts);

  ///
  /// exchange integer quantities
  ///
  idvectType idRecvBuff(newNoParts);
  idSendBuff.resize(noBuffParts);

  idvectPtrSetType::const_iterator intsItr = _quantities.ints.begin();
  idvectPtrSetType::const_iterator intsEnd = _quantities.ints.end();
  while (intsItr != intsEnd)
    {
      ///
      /// copy particles staying on the node into the buffer
      ///
      idvectRefType curInts(**intsItr);
      size_t storeIndex = stayOffset;
      stayItr = _partsIndices[myDomain].begin();
      while (stayItr != stayEnd)
        {
          idRecvBuff[storeIndex] = curInts[*stayItr];
          storeIndex++;
          stayItr++;
        }


      ///
      /// exchange non-staying particles
      ///
      queuedExch(**intsItr, idRecvBuff, partsQueue, localOffsets, noRecvParts);

      ///
      /// resize the quantitiy and copy buffer to particle manager
      ///
      PartManager.resize(**intsItr);
      idvectRangeType idLocal(**intsItr, newLocRange);
      idLocal = idRecvBuff;

      intsItr++;
    }

  ///
  /// exchange scalar quantities
  ///
  valvectType scalRecvBuff(newNoParts);
  scalSendBuff.resize(noBuffParts);

  valvectPtrSetType::const_iterator scalItr = _quantities.scalars.begin();
  valvectPtrSetType::const_iterator scalEnd = _quantities.scalars.end();
  while (scalItr != scalEnd)
    {
      ///
      /// copy particles staying on the node into the buffer
      ///
      valvectRefType curScal(**scalItr);
      size_t storeIndex = stayOffset;
      stayItr = _partsIndices[myDomain].begin();
      while (stayItr != stayEnd)
        {
          scalRecvBuff[storeIndex] = curScal[*stayItr];
          storeIndex++;
          stayItr++;
        }

      ///
      /// exchange non-staying particles
      ///
      queuedExch(**scalItr, scalRecvBuff, partsQueue,
                 localOffsets, noRecvParts);

      ///
      /// resize the quantitiy and copy buffer to particle manager
      ///
      PartManager.resize(**scalItr);
      valvectRangeType scalLocal(**scalItr, newLocRange);
      scalLocal = scalRecvBuff;

      scalItr++;
    }

  ///
  /// exchange vectorial quantities
  ///
  matrixType vectRecvBuff(newNoParts, 0);
  vectSendBuff.resize(noBuffParts, 0);

  matrixPtrSetType::const_iterator vectItr = _quantities.vects.begin();
  matrixPtrSetType::const_iterator vectEnd = _quantities.vects.end();
  while (vectItr != vectEnd)
    {
      ///
      /// copy particles staying on the node into the buffer
      ///
      matrixRefType curVect(**vectItr);

      ///
      /// resize 2nd dimension of buffers if necessary
      ///
      if (curVect.size2() != vectRecvBuff.size2())
        {
          vectRecvBuff.resize(newNoParts, curVect.size2(), false);
          vectSendBuff.resize(noBuffParts, curVect.size2(), false);
        }

      size_t storeIndex = stayOffset;
      stayItr = _partsIndices[myDomain].begin();
      while (stayItr != stayEnd)
        {
          particleRowType(vectRecvBuff, storeIndex) =
            particleRowType(curVect, *stayItr);
          storeIndex++;
          stayItr++;
        }

      ///
      /// exchange non-staying particles
      ///
      queuedExch(**vectItr, vectRecvBuff, partsQueue,
                 localOffsets, noRecvParts);

      ///
      /// resize the quantitiy and copy buffer to particle manager
      ///
      PartManager.resize(**vectItr);
      rangeType secDimRange(0, (**vectItr).size2());
      matrixRangeType vectLocal(**vectItr, newLocRange, secDimRange);
      vectLocal = vectRecvBuff;

      vectItr++;
    }

  // resize the rest of the variables
  PartManager.resizeAll();
}

void CommunicationManager::sendGhostsPrepare(
  domainPartsIndexRefType _ghostsIndices)
{
  prepareQueuesOffsets(_ghostsIndices, newNoParts, ghostQueue,
                       ghostOffsets, noRecvGhosts);
}

void CommunicationManager::sendGhosts(idvectRefType _idVect)
{
  queuedExch(_idVect, _idVect, ghostQueue, ghostOffsets, noRecvGhosts);
}

void CommunicationManager::sendGhosts(valvectRefType _valVect)
{
  queuedExch(_valVect, _valVect, ghostQueue, ghostOffsets, noRecvGhosts);
}

void CommunicationManager::sendGhosts(matrixRefType _matrix)
{
  queuedExch(_matrix, _matrix, ghostQueue, ghostOffsets, noRecvGhosts);
}

void CommunicationManager::sendGhosts(quantsTypeRefType _quantities)
{
  matrixPtrSetType::const_iterator vectItr = _quantities.vects.begin();
  matrixPtrSetType::const_iterator vectEnd = _quantities.vects.end();

  while (vectItr != vectEnd)
    {
      queuedExch(**vectItr, **vectItr,
                 ghostQueue, ghostOffsets, noRecvGhosts);
      vectItr++;
    }

  idvectPtrSetType::const_iterator intsItr = _quantities.ints.begin();
  idvectPtrSetType::const_iterator intsEnd = _quantities.ints.end();
  while (intsItr != intsEnd)
    {
      queuedExch(**intsItr, **intsItr,
                 ghostQueue, ghostOffsets, noRecvGhosts);
      intsItr++;
    }


  valvectPtrSetType::const_iterator scalItr = _quantities.scalars.begin();
  valvectPtrSetType::const_iterator scalEnd = _quantities.scalars.end();
  while (scalItr != scalEnd)
    {
      queuedExch(**scalItr, **scalItr,
                 ghostQueue, ghostOffsets, noRecvGhosts);
      scalItr++;
    }
}

///
/// generic function for chunked and queued exchange
///
template<class T>
void CommunicationManager::queuedExch(T& _src, T& _target,
                                      const sendQueueRefType _queue,
                                      countsVectType _offsets,
                                      countsVectType _noRecvParts)
{
  ///
  /// setup send queue iterators and send buffer
  ///
  sendQueueType::const_iterator qItr = _queue.begin();
  sendQueueType::const_iterator qEnd = _queue.end();

  ///
  /// setup receivers if anything's expected
  ///
  for (size_t i = 0; i < noDomains; i++)
    {
      if (_noRecvParts[i] > 0)
        {
          const size_t noEffRecvParts =
            std::min(noBuffParts, static_cast<size_t>(_noRecvParts[i]));
          recvChunk(_target, _offsets[i], recvReqs[i], noEffRecvParts, i);
        }
    }


  bool commFinished = false;
  while (!commFinished)
    {
      // Assume we have finished, prove otherwise
      commFinished = true;

      // check for arrived particles
      for (size_t i = 0; i < noDomains; i++)
        {
          if (_noRecvParts[i] > 0)
            {
              commFinished = false;
              if (recvReqs[i].Test())
                {
                  _noRecvParts[i] -= std::min(noBuffParts,
                                              static_cast<size_t>(_noRecvParts[i]));
                  // wire new receiver if we are still expecting particles
                  if (_noRecvParts[i] > 0)
                    {
                      _offsets[i] += noBuffParts;
                      const size_t noEffRecvParts =
                        std::min(noBuffParts,
                                 static_cast<size_t>(_noRecvParts[i]));
                      recvChunk(_target, _offsets[i],
                                recvReqs[i], noEffRecvParts, i);
                    }
                }
            }
        }

      ///
      /// send a send queue element
      ///
      if (qItr != qEnd)
        {
          commFinished = false;
          sendChunk(_src, *qItr);
          qItr++;
        }
    }
};

void CommunicationManager::sendChunk(const idvectRefType _src,
                                     sendQueueElemType _qelem)
{
  ///
  /// fill sendBuff
  ///
  static size_t count;

  count = 0;
  while (_qelem.itr != _qelem.end)
    {
      idSendBuff(count) = _src(*(_qelem.itr));
      _qelem.itr++;
      count++;
    }

  ///
  /// send chunk (may the return after a non-completed send be a problem,
  ///             when another send follows immediately trying to use
  ///             the same buffer?)
  ///
  MPI::COMM_WORLD.Ssend(&idSendBuff(0), count * identTypeSize, MPI::BYTE,
                        _qelem.sendRank, myRank);
};

void CommunicationManager::sendChunk(const valvectRefType _src,
                                     sendQueueElemType _qelem)
{
  ///
  /// fill sendBuff
  ///
  static size_t count;

  count = 0;
  while (_qelem.itr != _qelem.end)
    {
      scalSendBuff(count) = _src(*(_qelem.itr));
      _qelem.itr++;
      count++;
    }

  MPI::COMM_WORLD.Ssend(&scalSendBuff(0), count * valueTypeSize, MPI::BYTE,
                        _qelem.sendRank, myRank);
};

void CommunicationManager::sendChunk(const matrixRefType _src,
                                     sendQueueElemType _qelem)
{
  ///
  /// fill sendBuff
  ///
  static size_t count;

  count = 0;
  while (_qelem.itr != _qelem.end)
    {
      particleRowType(vectSendBuff, count) =
        particleRowType(_src, *(_qelem.itr));
      _qelem.itr++;
      count++;
    }

  MPI::COMM_WORLD.Ssend(&vectSendBuff(0, 0),
                        count * vectSendBuff.size2() * valueTypeSize,
                        MPI::BYTE,
                        _qelem.sendRank, myRank);
};


void CommunicationManager::recvChunk(idvectRefType _target,
                                     countsRefType _offset,
                                     MPI::Request& _recvRq,
                                     const size_t& _noParts,
                                     const size_t& _recvFrom)
{
  _recvRq = MPI::COMM_WORLD.Irecv(&_target(_offset), _noParts * identTypeSize,
                                  MPI::BYTE, _recvFrom, _recvFrom);
};

void CommunicationManager::recvChunk(valvectRefType _target,
                                     countsRefType _offset,
                                     MPI::Request& _recvRq,
                                     const size_t& _noParts,
                                     const size_t& _recvFrom)
{
  _recvRq = MPI::COMM_WORLD.Irecv(&_target(_offset), _noParts * valueTypeSize,
                                  MPI::BYTE, _recvFrom, _recvFrom);
};

void CommunicationManager::recvChunk(matrixRefType _target,
                                     countsRefType _offset,
                                     MPI::Request& _recvRq,
                                     const size_t& _noParts,
                                     const size_t& _recvFrom)
{
  _recvRq = MPI::COMM_WORLD.Irecv(&_target(_offset, 0), _noParts * valueTypeSize * _target.size2(),
                                  MPI::BYTE, _recvFrom, _recvFrom);
};

///
/// get from a particle index vector and a total offset constant the resulting
/// send queue, reveive offsets and number of receiving particles
///
void CommunicationManager::prepareQueuesOffsets(domainPartsIndexRefType _indices,
                                                size_t _totalOffset,
                                                sendQueueRefType _queue,
                                                countsVectRefType _offsets,
                                                countsVectRefType _noRecvParts)
{
  const size_t noRanks = noDomains;

  ///
  /// determine _noRecvParts
  ///
  countsVectType noSendParts(noDomains);

  for (size_t i = 0; i < noDomains; i++)
    {
      noSendParts[ domainToRank[i] ] = _indices[i].size();
    }
  MPI::COMM_WORLD.Alltoall(&noSendParts[0], 1, MPI::INT,
                           &_noRecvParts[0], 1, MPI::INT);

  ///
  /// determine offsets
  ///
  _offsets[0] = _totalOffset;
  for (size_t i = 1; i < noDomains; i++)
    {
      _offsets[i] = _offsets[i - 1] + _noRecvParts[i - 1];
    }

  ///
  /// fill queue
  ///
  _queue.clear();
  for (size_t i = 1; i < noRanks; i++)
    {
      const size_t curSendRank = (myRank + i) % noRanks;

      partsIndexVectType::const_iterator idxItr =
        _indices[ rankToDomain[ curSendRank ] ].begin();
      partsIndexVectType::const_iterator lastIndex =
        _indices[ rankToDomain[ curSendRank ] ].end();

      while (idxItr != lastIndex)
        {
          sendQueueElemType newElem;
          newElem.sendRank = curSendRank;
          newElem.itr = idxItr;

          static size_t noQParts;
          noQParts = 0;
          while (idxItr != lastIndex && noQParts < noBuffParts)
            {
              idxItr++;
              noQParts++;
            }
          newElem.end = idxItr;
          newElem.noParts = noQParts;

          _queue.push_back(newElem);
        }
    }
}

void CommunicationManager::sendMatrix(matrixRefType _matrix, size_t _recvDomain)
{
  static size_t noRemElems, noCurElems, noCurBytes, round;
  const size_t noBuffElems = commBuffSize / sizeof(valueType);

  const size_t recvRank = domainToRank[_recvDomain];

  ///
  /// no checks are performed whether the matrix
  /// on the receiving side actually matches the amount
  /// of incoming data!
  ///
  noRemElems = _matrix.size1() * _matrix.size2();
  round = 0;

  while (noRemElems > 0)
    {
      noCurElems = std::min(noRemElems, noBuffElems);
      noCurBytes = noCurElems * sizeof(valueType);

      ///
      /// this pointer arithmetics stuff assumes continous
      /// storage for matrixType
      ///
      MPI::COMM_WORLD.Send((&_matrix(0, 0) + round * noBuffElems),
                           noCurBytes, MPI_BYTE, recvRank, myRank + round);
      noRemElems -= noCurElems;
      round++;
    }
}

void CommunicationManager::recvMatrix(matrixRefType _matrix, size_t _sendDomain)
{
  static size_t noRemElems, noCurElems, noCurBytes, round;
  const size_t noBuffElems = commBuffSize / sizeof(valueType);

  const size_t sendRank = domainToRank[_sendDomain];

  ///
  /// no checks are performed whether the matrix
  /// on the receiving side actually matches the amount
  /// of incoming data!
  ///
  noRemElems = _matrix.size1() * _matrix.size2();
  round = 0;

  while (noRemElems > 0)
    {
      noCurElems = std::min(noRemElems, noBuffElems);
      noCurBytes = noCurElems * sizeof(valueType);

      ///
      /// this pointer arithmetics stuff assumes continous
      /// storage for matrixType
      ///
      MPI::COMM_WORLD.Recv((&_matrix(0, 0) + round * noBuffElems),
                           noCurBytes, MPI_BYTE, sendRank, sendRank + round);
      noRemElems -= noCurElems;
      round++;
    }
}

void CommunicationManager::sendBitset(bitsetRefType _bitset, size_t _recvRank)
{
  /// noBlock = no. of bits / no. of bits per bitsetBlockType
  size_t noBlocks = lrint(
    ceil(static_cast<double>(_bitset.size())
         / static_cast<double>(sizeof(bitsetBlockType) * 8.)));

  ///
  /// allocate buffer and fill it
  ///
  std::vector<bitsetBlockType> sendBuff(noBlocks);
  boost::to_block_range(_bitset, sendBuff.begin());

  ///
  /// send it
  ///
  sendVector<bitsetBlockType>(sendBuff, _recvRank);
}

void CommunicationManager::recvBitset(bitsetRefType _bitset, size_t _sendDomain)
{
  /// noBlock = no. of bits / no. of bits per bitsetBlockType
  size_t noBlocks = lrint(
    ceil(static_cast<double>(_bitset.size())
         / static_cast<double>(sizeof(bitsetBlockType) * 8.)));

  ///
  /// allocate buffer
  ///
  std::vector<bitsetBlockType> recvBuff(noBlocks);

  ///
  /// send it
  ///
  recvVector<bitsetBlockType>(recvBuff, _sendDomain);

  ///
  /// copy buffer to bitset
  ///
  boost::from_block_range(recvBuff.begin(), recvBuff.end(), _bitset);
}

void CommunicationManager::sumUpCounts(countsVectRefType _indexVect)
{
  static size_t noRemElems, noCurElems, round;
  const size_t noBuffElems = commBuffSize / sizeof(countsType);

  static countsVectType rootBuff;

  noRemElems = _indexVect.size();
  round = 0;

  while (noRemElems > 0)
    {
      noCurElems = std::min(noRemElems, noBuffElems);
      ///
      /// check whether this MPI implementation supports
      /// the faster MPI_IN_PLACE method
      ///
#ifdef MPI_IN_PLACE
      ///
      /// the vector cast is ugly, but actually legal:
      /// http://www.parashift.com/c++-faq-lite/containers.html#faq-34.3
      ///
      MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,
                                &_indexVect[round * noBuffElems],
                                noCurElems, MPI::INT, MPI::SUM);
#else
      ///
      /// slower without MPI_IN_PLACE
      ///
      static countsVectType recvBuffer;
      recvBuffer.resize(noCurElems);
      MPI::COMM_WORLD.Allreduce(&_indexVect[round * noBuffElems],
                                &recvBuffer[0],
                                noCurElems, MPI::INT, MPI::SUM);
      const size_t offset = round * noBuffElems;
      for (size_t i = 0; i < noCurElems; i++)
        {
          _indexVect[i + offset] = recvBuffer[i];
        }
#endif
      noRemElems -= noCurElems;
      round++;
    }
}

void CommunicationManager::max(valueRefType _val)
{
  static double doubleBuff;

  doubleBuff = static_cast<double>(_val);

#ifdef MPI_IN_PLACE
  MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, &doubleBuff, 1,
                            MPI::DOUBLE, MPI::MAX);
#else
  static double recvDoubleBuff;
  MPI::COMM_WORLD.Allreduce(&doubleBuff, &recvDoubleBuff, 1,
                            MPI::DOUBLE, MPI::MAX);
  doubleBuff = recvDoubleBuff;
#endif
  _val = static_cast<valueType>(doubleBuff);
}

void CommunicationManager::min(valueRefType _val)
{
  static double doubleBuff;

  doubleBuff = static_cast<double>(_val);

#ifdef MPI_IN_PLACE
  MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, &doubleBuff, 1,
                            MPI::DOUBLE, MPI::MIN);
#else
  static double recvDoubleBuff;
  MPI::COMM_WORLD.Allreduce(&doubleBuff, &recvDoubleBuff, 1,
                            MPI::DOUBLE, MPI::MIN);
  doubleBuff = recvDoubleBuff;
#endif
  _val = static_cast<valueType>(doubleBuff);
}

void CommunicationManager::sum(valueRefType _val)
{
  static double doubleBuff;

  doubleBuff = static_cast<double>(_val);

#ifdef MPI_IN_PLACE
  MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, &doubleBuff, 1,
                            MPI::DOUBLE, MPI::SUM);
#else
  static double recvDoubleBuff;
  MPI::COMM_WORLD.Allreduce(&doubleBuff, &recvDoubleBuff, 1,
                            MPI::DOUBLE, MPI::SUM);
  doubleBuff = recvDoubleBuff;
#endif
  _val = static_cast<valueType>(doubleBuff);
}

template<class T>
void CommunicationManager::sendVector(std::vector<T>& _vector,
                                      size_t _recvDomain)
{
  static size_t noRemElems, noCurElems, noCurBytes, round;
  const size_t noBuffElems = commBuffSize / sizeof(T);

  const size_t recvRank = domainToRank[_recvDomain];

  ///
  /// no checks are performed whether the vector
  /// on the receiving side actually matches the amount
  /// of incoming data!
  ///
  noRemElems = _vector.size();
  round = 0;

  while (noRemElems > 0)
    {
      noCurElems = std::min(noRemElems, noBuffElems);
      noCurBytes = noCurElems * sizeof(T);

      ///
      /// this pointer arithmetics stuff assumes continous
      /// storage for std::vector<T>
      ///
      /// does not work for std::vector<bool> !
      ///
      MPI::COMM_WORLD.Send((&_vector[0] + round * noBuffElems),
                           noCurBytes, MPI_BYTE, recvRank,
                           myRank + 255 + round);
      noRemElems -= noCurElems;
      round++;
    }
}

template<class T>
void CommunicationManager::recvVector(std::vector<T>& _vector,
                                      size_t _sendDomain)
{
  static size_t noRemElems, noCurElems, noCurBytes, round;
  const size_t noBuffElems = commBuffSize / sizeof(T);

  const size_t sendRank = domainToRank[_sendDomain];

  ///
  /// no checks are performed whether the vector
  /// on the receiving side actually matches the amount
  /// of incoming data!
  ///
  noRemElems = _vector.size();
  round = 0;

  while (noRemElems > 0)
    {
      noCurElems = std::min(noRemElems, noBuffElems);
      noCurBytes = noCurElems * sizeof(T);

      ///
      /// this pointer arithmetics stuff assumes continous
      /// storage for std::vector<T>
      ///
      /// does not work for std::vector<bool> !
      ///
      MPI::COMM_WORLD.Recv((&_vector[0] + round * noBuffElems),
                           noCurBytes, MPI_BYTE, sendRank,
                           sendRank + 255 + round);
      noRemElems -= noCurElems;
      round++;
    }
}

size_t CommunicationManager::getMyDomain()
{
  return myDomain;
}
size_t CommunicationManager::getNoDomains()
{
  return noDomains;
}

void CommunicationManager::barrier()
{
  MPI::COMM_WORLD.Barrier();
}

size_t CommunicationManager::domainToMPIrank(size_t _domain)
{
  return domainToRank[_domain];
}
size_t CommunicationManager::MPIrankToDomain(size_t _rank)
{
  return rankToDomain[_rank];
}
};

#endif
