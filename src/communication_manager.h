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

void exchange(domainPartsIndexRefType _partsIndices,
              domainPartsIndexRefType _ghostsIndices,
              quantsTypeRefType _quantities);

void sendGhosts(idvectRefType _idVect);
void sendGhosts(valvectRefType _valVect);
void sendGhosts(matrixRefType _matrix);

private:
std::vector<MPI::Request> recvReqs;
template<class T> void queuedExch(T& _src, T& _trgt,
                                  const sendQueueRefType _queue,
                                  countsVectType _offsets,
                                  countsVectType _noRecvParts);

void sendChunk(const idvectRefType _src, sendQueueElemType _qelem);
idvectType idSendBuff;

//void sendChunk(const valvectRefType _src, sendQueueElemType _qelem);

void recvChunk(idvectRefType _target, countsRefType _offset, MPI::Request& _recvRq,
               const size_t& _noParts, const size_t& _recvFrom);
/*void recvChunk(valvectRefType _target, size_t& _offset, MPI::Request& _recvRq,
               size_t& _noParts, size_t& _recvFrom);*/

void prepareQueuesOffsets(domainPartsIndexRefType _indices, size_t _totalOffset,
                          sendQueueRefType _queue,
                          countsVectRefType _offsets,
                          countsVectRefType _noRecvParts);

sendQueueType ghostQueue;
countsVectType localOffsets, noRecvParts, ghostOffsets, noRecvGhosts;
size_t identTypeSize, valueTypeSize;

public:
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

std::vector<matrixType> composeBuffers;
std::vector<matrixType> recvBuffers;
matrixType sendBuffer;

std::vector<size_t> sendNoPart;
std::vector<size_t> recvNoPart;

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
//size_t CommunicationManager::commBuffSize = 131072; //  128 kByte
size_t CommunicationManager::commBuffSize = 1048576; // 1024 kByte

///
/// buffer size in particles for the exchange() functions
/// 64k gives about 1.5MB for a 3D double precision vector
///
/// you may want to lower this value, when you plan to transmit higher dimension
/// vectors like tensors
///
size_t CommunicationManager::noBuffParts = 65536;
//size_t CommunicationManager::noBuffParts = 200; // for testing purposes

void CommunicationManager::exchange(domainPartsIndexRefType _partsIndices,
                                    domainPartsIndexRefType _ghostsIndices,
                                    quantsTypeRefType _quantities)
{
  // exchange number of particles

  double startTime = MPI_Wtime();
  sendQueueType partsQueue;

  prepareQueuesOffsets(_partsIndices, 0, partsQueue,
                       localOffsets, noRecvParts);

  ///
  /// prepare particle queue and number of particles
  ///
  size_t noParts = 0;
  for (size_t i = 0; i < noDomains; i++)
    {
      noParts += noRecvParts[i];
    }

  rangeType newLocRange(0, noParts);

  std::cout << myDomain << ": " << std::fixed << std::right << std::setw(15) << std::setprecision(6)
            << MPI_Wtime() - startTime
            << " prepared partsQueue " << "\n";

  ///
  /// prepare ghost queue and number of ghosts
  ///
  prepareQueuesOffsets(_ghostsIndices, noParts, ghostQueue,
                       ghostOffsets, noRecvGhosts);

  size_t noGhosts = 0;
  for (size_t i = 0; i < noDomains; i++)
    {
      noGhosts += noRecvGhosts[i];
    }

  std::cout << myDomain << ": " << std::fixed << std::right << std::setw(15) << std::setprecision(6)
            << MPI_Wtime() - startTime
            << " prepared ghostQueue " << "\n";

  ///
  /// determine information for particles staying on node
  ///
  const size_t stayOffset = localOffsets[myRank];
  noRecvParts[myRank] = 0;

  partsIndexVectType::const_iterator stayItr;
  partsIndexVectType::const_iterator stayEnd = _partsIndices[myDomain].end();

  ///
  /// set new number of particles
  ///
  PartManager.setNoParts(noParts, noGhosts);

  /// exchange ints

  // set up buffers
  idvectType idRecvBuff(noParts);
  idSendBuff.resize(noBuffParts);

  std::cout << myDomain << ": " << std::fixed << std::right << std::setw(15) << std::setprecision(6)
            << MPI_Wtime() - startTime
            << " set up buffers      " << "\n";

  idvectPtrSetType::const_iterator intsItr = _quantities.ints.begin();
  idvectPtrSetType::const_iterator intsEnd = _quantities.ints.end();
  //exchangeInt();
  while ( intsItr != intsEnd )
  {
    //idvectRefType id(**intsItr);

    ///
    /// move local data to its designated space
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
    //std::cout << myDomain << ": " << stayOffset << " - " << storeIndex << "    local\n";

    std::cout << myDomain << ": " << std::fixed << std::right << std::setw(15) << std::setprecision(6)
              << MPI_Wtime() - startTime
              << " copied local stuff  " << "\n";

    queuedExch(**intsItr, idRecvBuff, partsQueue, localOffsets, noRecvParts);

    std::cout << myDomain << ": " << std::fixed << std::right << std::setw(15) << std::setprecision(6)
              << MPI_Wtime() - startTime
              << " exchanged remote stuff " << "\n";

    //std::cout << myDomain << ": size of new id is " << idRecvBuff.size() << "\n";

    // set no. of particles
    PartManager.resize(**intsItr);
    //std::cout << myDomain << ": " << id.size() << "\n";

    std::cout << myDomain << ": " << std::fixed << std::right << std::setw(15) << std::setprecision(6)
              << MPI_Wtime() - startTime
              << " resized container   " << "\n";

    idvectRangeType idLocal(**intsItr, newLocRange);
    idLocal = idRecvBuff;

    std::cout << myDomain << ": " << std::fixed << std::right << std::setw(15) << std::setprecision(6)
              << MPI_Wtime() - startTime
              << " copied back local buff\n" << "\n";
    
    intsItr++;
    //std::cout << myDomain << ": " << idLocal.size() << "\n";
  }

  // exchange scalars
  // set up buffers
  //exchangeScalar();
  // set no. of particles
  // copy buffers

  // exchange vectors
  // set up buffers
  //exchangeVect();
  // set no. of particles
  // copy buffers

  // resize the rest of the variables
  //PartManager.resizeAll();
}

void CommunicationManager::sendGhosts(idvectRefType _idVect)
{
  queuedExch(_idVect, _idVect, ghostQueue, ghostOffsets, noRecvGhosts);
}

void CommunicationManager::sendGhosts(valvectRefType _valVect)
{
  //exchangeScalar(_valVect, _valVect, ghostQueue, ghostOffsets, noRecvGhosts);
}

void CommunicationManager::sendGhosts(matrixRefType _matrix)
{
  //exchangeVect(_matrix, _matrix, ghostQueue, ghostOffsets, noRecvGhosts);
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
      //if (myDomain == 0)
      //  std::cout << count << " " << *(_qelem.itr) << "\n";
      idSendBuff(count) = _src(*(_qelem.itr));
      /// somethings fishy: on node33 (gcc 4.1) the following line just increments the iterators value, but not the iterator itself!
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

void CommunicationManager::recvChunk(idvectRefType _target,
                                     countsRefType _offset,
                                     MPI::Request& _recvRq,
                                     const size_t& _noParts,
                                     const size_t& _recvFrom)
{
  //std::cout << myDomain << ": " << _offset << " - " << _offset+_noParts << "    from rank " << _recvFrom << "\n";
  _recvRq = MPI::COMM_WORLD.Irecv(&_target(_offset), _noParts * identTypeSize,
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
  MPI::COMM_WORLD.Alltoall(&noSendParts[0], 1, MPI::UNSIGNED_LONG,
                           &_noRecvParts[0], 1, MPI::UNSIGNED_LONG);

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

          static size_t noParts;
          noParts = 0;
          while (idxItr != lastIndex && noParts < noBuffParts)
            {
              idxItr++;
              noParts++;
            }
          newElem.end = idxItr;
          newElem.noParts = noParts;

          _queue.push_back(newElem);
        }
    }
}

//sendQueueType      ghostQueue;
//countsVectType localOffsets, ghostOffsets;

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
  MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, &doubleBuff, 1, MPI::DOUBLE, MPI::MAX);
#else
  static double recvDoubleBuff;
  MPI::COMM_WORLD.Allreduce(&doubleBuff, &recvDoubleBuff, 1, MPI::DOUBLE, MPI::MAX);
  doubleBuff = recvDoubleBuff;
#endif
  _val = static_cast<valueType>(doubleBuff);
}

void CommunicationManager::min(valueRefType _val)
{
  static double doubleBuff;

  doubleBuff = static_cast<double>(_val);

#ifdef MPI_IN_PLACE
  MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, &doubleBuff, 1, MPI::DOUBLE, MPI::MIN);
#else
  static double recvDoubleBuff;
  MPI::COMM_WORLD.Allreduce(&doubleBuff, &recvDoubleBuff, 1, MPI::DOUBLE, MPI::MIN);
  doubleBuff = recvDoubleBuff;
#endif
  _val = static_cast<valueType>(doubleBuff);
}

void CommunicationManager::sum(valueRefType _val)
{
  static double doubleBuff;

  doubleBuff = static_cast<double>(_val);

#ifdef MPI_IN_PLACE
  MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, &doubleBuff, 1, MPI::DOUBLE, MPI::SUM);
#else
  static double recvDoubleBuff;
  MPI::COMM_WORLD.Allreduce(&doubleBuff, &recvDoubleBuff, 1, MPI::DOUBLE, MPI::SUM);
  doubleBuff = recvDoubleBuff;
#endif
  _val = static_cast<valueType>(doubleBuff);
}

template<class T>
void CommunicationManager::sendVector(std::vector<T>& _vector, size_t _recvDomain)
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
                           noCurBytes, MPI_BYTE, recvRank, myRank + 255 + round);
      noRemElems -= noCurElems;
      round++;
    }
}

template<class T>
void CommunicationManager::recvVector(std::vector<T>& _vector, size_t _sendDomain)
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
                           noCurBytes, MPI_BYTE, sendRank, sendRank + 255 + round);
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
