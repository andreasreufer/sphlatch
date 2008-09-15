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

typedef sphlatch::ParticleManager PartManagerType;

static self_reference instance(void);

PartManagerType& PartManager;

///
/// this is the central exchange part
///
void exchange(domainPartsIndexRefType _partsIndices, size_t _noGhosts);

quantsType exchangeQuants;

void sendGhostsPrepare(domainPartsIndexRefType _ghostsIndices);

void sendGhosts(idvectRefType _idVect);
void sendGhosts(valvectRefType _valVect);
void sendGhosts(matrixRefType _matrix);

void sendGhosts(quantsRefType _quantities);

void regExchQuant(idvectRefType _idVect);
void regExchQuant(valvectRefType _valVect);
void regExchQuant(matrixRefType _matrix);

private:
void calcOffsets(domainPartsIndexRefType _indices,
                 size_t _totalOffset,
                 countsVectRefType _sendOffset,
                 countsVectRefType _partsTo,
                 countsVectRefType _recvOffset,
                 countsVectRefType _partsFrom);

domainPartsIndexType ghostIndices;
countsVectType gSendOffsets, gRecvOffsets, gPartsFrom, gPartsTo;
size_t identSize, floatSize, newNoParts, noSendGhosts;

///
/// various little helper functions
///
///
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
double wtime();

size_t domainToMPIrank(size_t _domain);
size_t MPIrankToDomain(size_t _rank);

protected:
CommunicationManager(void);
~CommunicationManager(void);

private:
static self_pointer _instance;
static size_t commBuffSize;
size_t myDomain, myRank, noDomains, noRanks;

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
  noRanks = noDomains;
  domainToRank.resize(noDomains);
  rankToDomain.resize(noDomains);

  ///
  /// insert mapping here
  ///
  /// right now rank maps directly to the domain number
  ///
  for (size_t i = 0; i < noDomains; i++)
    {
      domainToRank[i] = i;
      rankToDomain[i] = i;
    }

  myRank = MPI::COMM_WORLD.Get_rank();
  myDomain = rankToDomain[ myRank ];

  gSendOffsets.resize(noDomains);
  gRecvOffsets.resize(noDomains);
  gPartsFrom.resize(noDomains);
  gPartsTo.resize(noDomains);

  identSize = sizeof(identType);
  floatSize = sizeof(valueType);
}

CommunicationManager::~CommunicationManager(void)
{
}

///
/// buffer size in bytes, MPI transfers will have on average this size
///
size_t CommunicationManager::commBuffSize = 131072;   //  128 kByte
//size_t CommunicationManager::commBuffSize = 1048576;    // 1024 kByte
//size_t CommunicationManager::commBuffSize = 8388608;  // 8096 kByte

void CommunicationManager::exchange(domainPartsIndexRefType _partsIndices,
                                    size_t _noGhosts)
{
  ///
  /// prepare particle queue and number of particles
  ///
  countsVectType sendOffset, recvOffset, partsTo, partsFrom;

  sendOffset.resize(noDomains);
  partsTo.resize(noDomains);
  recvOffset.resize(noDomains);
  partsFrom.resize(noDomains);

  calcOffsets(_partsIndices, 0, sendOffset, partsTo, recvOffset, partsFrom);

  newNoParts = 0;
  size_t noSendParts = 0;
  for (size_t i = 0; i < noRanks; i++)
    {
      newNoParts += partsFrom[i];
      noSendParts += partsTo[i];
    }

  ///
  /// current number of particles
  ///
  //const valueType noParts = PartManager.getNoLocalParts();

  ///
  /// set new number of particles (noGhosts needs to be known)
  ///
  PartManager.setNoParts(newNoParts, _noGhosts);

  ///
  /// iterators for the particle index vectors
  ///
  partsIndexVectType::const_iterator idxItr, idxEnd;

  ///
  /// vectors for byte counts and offsets
  ///
  countsVectType sendByteOffset, recvByteOffset, bytesTo, bytesFrom;
  sendByteOffset.resize(noDomains);
  recvByteOffset.resize(noDomains);
  bytesTo.resize(noDomains);
  bytesFrom.resize(noDomains);

  ///
  /// exchange integer quantities
  ///
  idvectType idSendBuff(noSendParts);

  idvectPtrSetType::const_iterator intsItr = exchangeQuants.ints.begin();
  idvectPtrSetType::const_iterator intsEnd = exchangeQuants.ints.end();
  while (intsItr != intsEnd)
    {
      idvectRefType curInts(**intsItr);

      ///
      /// fill send buffer and prepare byte count vectors
      ///
      for (size_t i = 0; i < noDomains; i++)
        {
          const size_t curRank = domainToRank[i];
          size_t buffIdx = sendOffset[curRank];

          idxItr = _partsIndices[i].begin();
          idxEnd = _partsIndices[i].end();

          while (idxItr != idxEnd)
            {
              idSendBuff(buffIdx) = curInts(*idxItr);
              buffIdx++;
              idxItr++;
            }

          sendByteOffset[i] = identSize * sendOffset[i];
          recvByteOffset[i] = identSize * recvOffset[i];
          bytesTo[i] = identSize * partsTo[i];
          bytesFrom[i] = identSize * partsFrom[i];
        }

      ///
      /// resize the quantitiy
      ///
      PartManager.resize(**intsItr, false);

      ///
      /// exchange data
      ///
      barrier();
      MPI::COMM_WORLD.Alltoallv(&idSendBuff(0), &bytesTo[0],
                                &sendByteOffset[0], MPI::BYTE,
                                &curInts(0), &bytesFrom[0],
                                &recvByteOffset[0], MPI::BYTE);
      intsItr++;
    }
  idSendBuff.resize(0, false);

  ///
  /// exchange scalar quantities
  ///
  valvectType scalSendBuff(noSendParts);

  valvectPtrSetType::const_iterator scalItr = exchangeQuants.scalars.begin();
  valvectPtrSetType::const_iterator scalEnd = exchangeQuants.scalars.end();
  while (scalItr != scalEnd)
    {
      valvectRefType curScal(**scalItr);

      ///
      /// fill send buffer and prepare byte count vectors
      ///
      for (size_t i = 0; i < noDomains; i++)
        {
          const size_t curRank = domainToRank[i];
          size_t buffIdx = sendOffset[curRank];

          idxItr = _partsIndices[i].begin();
          idxEnd = _partsIndices[i].end();

          while (idxItr != idxEnd)
            {
              scalSendBuff(buffIdx) = curScal(*idxItr);
              buffIdx++;
              idxItr++;
            }

          sendByteOffset[i] = floatSize * sendOffset[i];
          recvByteOffset[i] = floatSize * recvOffset[i];
          bytesTo[i] = floatSize * partsTo[i];
          bytesFrom[i] = floatSize * partsFrom[i];
        }

      ///
      /// resize the quantitiy
      ///
      PartManager.resize(**scalItr, false);

      ///
      /// exchange data
      ///
      barrier();
      MPI::COMM_WORLD.Alltoallv(&scalSendBuff(0), &bytesTo[0],
                                &sendByteOffset[0], MPI::BYTE,
                                &curScal(0), &bytesFrom[0],
                                &recvByteOffset[0], MPI::BYTE);
      scalItr++;
    }
  scalSendBuff.resize(0, false);

  ///
  /// exchange vectorial quantities
  ///
  matrixType vectSendBuff(noSendParts, 0);

  matrixPtrSetType::const_iterator vectItr = exchangeQuants.vects.begin();
  matrixPtrSetType::const_iterator vectEnd = exchangeQuants.vects.end();
  while (vectItr != vectEnd)
    {
      matrixRefType curVect(**vectItr);
      const size_t matrSize2 = curVect.size2();

      ///
      /// resize 2nd dimension of buffers if necessary
      ///
      if (matrSize2 != vectSendBuff.size2())
        vectSendBuff.resize(vectSendBuff.size1(), matrSize2, false);

      ///
      /// fill send buffer and prepare byte count vectors
      ///
      for (size_t i = 0; i < noDomains; i++)
        {
          const size_t curRank = domainToRank[i];
          size_t buffIdx = sendOffset[curRank];

          idxItr = _partsIndices[i].begin();
          idxEnd = _partsIndices[i].end();

          while (idxItr != idxEnd)
            {
              particleRowType(vectSendBuff, buffIdx) =
                particleRowType(curVect, *idxItr);
              buffIdx++;
              idxItr++;
            }

          sendByteOffset[i] = floatSize * matrSize2 * sendOffset[i];
          recvByteOffset[i] = floatSize * matrSize2 * recvOffset[i];
          bytesTo[i] = floatSize * matrSize2 * partsTo[i];
          bytesFrom[i] = floatSize * matrSize2 * partsFrom[i];
        }

      ///
      /// resize the quantitiy
      ///
      PartManager.resize(**vectItr, false);

      ///
      /// exchange data
      ///
      barrier();
      MPI::COMM_WORLD.Alltoallv(&vectSendBuff(0, 0), &bytesTo[0],
                                &sendByteOffset[0], MPI::BYTE,
                                &curVect(0, 0), &bytesFrom[0],
                                &recvByteOffset[0], MPI::BYTE);
      vectItr++;
    }
  vectSendBuff.resize(0, 0, false);

  // resize the rest of the variables
  PartManager.resizeAll();
}

void CommunicationManager::sendGhostsPrepare(
  domainPartsIndexRefType _ghostsIndices)
{
  ghostIndices = _ghostsIndices;
  calcOffsets(ghostIndices, newNoParts,
              gSendOffsets, gPartsTo,
              gRecvOffsets, gPartsFrom);

  noSendGhosts = 0;
  for (size_t i = 0; i < noRanks; i++)
    {
      noSendGhosts += gPartsTo[i];
    }
}

void CommunicationManager::sendGhosts(idvectRefType _idVect)
{
  idvectType idSendBuff(noSendGhosts);

  countsVectType sendByteOffset(noDomains);
  countsVectType recvByteOffset(noDomains);
  countsVectType bytesTo(noDomains);
  countsVectType bytesFrom(noDomains);
  
  ///
  /// fill send buffer and prepare byte count vectors
  ///
  for (size_t i = 0; i < noDomains; i++)
    {
      const size_t curRank = domainToRank[i];
      size_t buffIdx = gSendOffsets[curRank];

      partsIndexVectType::const_iterator idxItr = ghostIndices[i].begin();
      partsIndexVectType::const_iterator idxEnd = ghostIndices[i].end();

      while (idxItr != idxEnd)
        {
          idSendBuff(buffIdx) = _idVect(*idxItr);
          buffIdx++;
          idxItr++;
        }

      sendByteOffset[i] = identSize * gSendOffsets[i];
      recvByteOffset[i] = identSize * gRecvOffsets[i];
      bytesTo[i] = identSize * gPartsTo[i];
      bytesFrom[i] = identSize * gPartsFrom[i];
    }

  ///
  /// exchange data
  ///
  barrier();
  MPI::COMM_WORLD.Alltoallv(&idSendBuff(0), &bytesTo[0],
                            &sendByteOffset[0], MPI::BYTE,
                            &_idVect(0), &bytesFrom[0],
                            &recvByteOffset[0], MPI::BYTE);
  barrier();
}

void CommunicationManager::sendGhosts(valvectRefType _valVect)
{
  valvectType valSendBuff(noSendGhosts);

  countsVectType sendByteOffset(noDomains);
  countsVectType recvByteOffset(noDomains);
  countsVectType bytesTo(noDomains);
  countsVectType bytesFrom(noDomains);
  
  ///
  /// fill send buffer and prepare byte count vectors
  ///
  for (size_t i = 0; i < noDomains; i++)
    {
      const size_t curRank = domainToRank[i];
      size_t buffIdx = gSendOffsets[curRank];

      partsIndexVectType::const_iterator idxItr = ghostIndices[i].begin();
      partsIndexVectType::const_iterator idxEnd = ghostIndices[i].end();

      while (idxItr != idxEnd)
        {
          valSendBuff(buffIdx) = _valVect(*idxItr);
          buffIdx++;
          idxItr++;
        }

      sendByteOffset[i] = floatSize * gSendOffsets[i];
      recvByteOffset[i] = floatSize * gRecvOffsets[i];
      bytesTo[i] = floatSize * gPartsTo[i];
      bytesFrom[i] = floatSize * gPartsFrom[i];
    }

  ///
  /// exchange data
  ///
  barrier();
  MPI::COMM_WORLD.Alltoallv(&valSendBuff(0), &bytesTo[0],
                            &sendByteOffset[0], MPI::BYTE,
                            &_valVect(0), &bytesFrom[0],
                            &recvByteOffset[0], MPI::BYTE);
  barrier();
}

void CommunicationManager::sendGhosts(matrixRefType _matrix)
{
  const size_t matrSize2 = _matrix.size2();
  matrixType vectSendBuff(noSendGhosts, matrSize2);
  
  countsVectType sendByteOffset(noDomains);
  countsVectType recvByteOffset(noDomains);
  countsVectType bytesTo(noDomains);
  countsVectType bytesFrom(noDomains);
  
  ///
  /// fill send buffer and prepare byte count vectors
  ///
  for (size_t i = 0; i < noDomains; i++)
    {
      const size_t curRank = domainToRank[i];
      size_t buffIdx = gSendOffsets[curRank];

      partsIndexVectType::const_iterator idxItr = ghostIndices[i].begin();
      partsIndexVectType::const_iterator idxEnd = ghostIndices[i].end();

      while (idxItr != idxEnd)
        {
          particleRowType(vectSendBuff, buffIdx) =
                particleRowType(_matrix, *idxItr);
          buffIdx++;
          idxItr++;
        }

      sendByteOffset[i] = floatSize * matrSize2 * gSendOffsets[i];
      recvByteOffset[i] = floatSize * matrSize2 * gRecvOffsets[i];
      bytesTo[i] = floatSize * matrSize2 * gPartsTo[i];
      bytesFrom[i] = floatSize * matrSize2 * gPartsFrom[i];
    }
  
  ///
  /// exchange data
  ///
  barrier();
  MPI::COMM_WORLD.Alltoallv(&vectSendBuff(0,0), &bytesTo[0],
                            &sendByteOffset[0], MPI::BYTE,
                            &_matrix(0,0), &bytesFrom[0],
                            &recvByteOffset[0], MPI::BYTE);
  barrier();
}

void CommunicationManager::sendGhosts(quantsRefType _quantities)
{
  matrixPtrSetType::const_iterator vectItr = _quantities.vects.begin();
  matrixPtrSetType::const_iterator vectEnd = _quantities.vects.end();

  while (vectItr != vectEnd)
    {
      sendGhosts(**vectItr);
      vectItr++;
    }

  idvectPtrSetType::const_iterator intsItr = _quantities.ints.begin();
  idvectPtrSetType::const_iterator intsEnd = _quantities.ints.end();
  while (intsItr != intsEnd)
    {
      sendGhosts(**intsItr);
      intsItr++;
    }


  valvectPtrSetType::const_iterator scalItr = _quantities.scalars.begin();
  valvectPtrSetType::const_iterator scalEnd = _quantities.scalars.end();
  while (scalItr != scalEnd)
    {
      sendGhosts(**scalItr);
      scalItr++;
    }
}

///
/// register a quantity for exchange (only when ptr is not already
/// registered)
///
void CommunicationManager::regExchQuant(idvectRefType _idVect)
{
  idvectPtrSetType::iterator searchItr = exchangeQuants.ints.find(&_idVect);

  if (searchItr == exchangeQuants.ints.end())
    {
      exchangeQuants.ints.insert(&_idVect);
    }
}

void CommunicationManager::regExchQuant(valvectRefType _valVect)
{
  valvectPtrSetType::iterator searchItr =
    exchangeQuants.scalars.find(&_valVect);

  if (searchItr == exchangeQuants.scalars.end())
    {
      exchangeQuants.scalars.insert(&_valVect);
    }
}

void CommunicationManager::regExchQuant(matrixRefType _matrix)
{
  matrixPtrSetType::iterator searchItr =
    exchangeQuants.vects.find(&_matrix);

  if (searchItr == exchangeQuants.vects.end())
    {
      exchangeQuants.vects.insert(&_matrix);
    }
}

///
/// calculate from a particle index vector and a total offset constant the resulting
/// counts and offsets for sending and receiving
///
/// MPI ranks are used insted of domain numbers for the indices of the offset
/// and count vectors!
///
void CommunicationManager::calcOffsets(domainPartsIndexRefType _indices,
                                       size_t _totalOffset,
                                       countsVectRefType _sendOffset,
                                       countsVectRefType _partsTo,
                                       countsVectRefType _recvOffset,
                                       countsVectRefType _partsFrom)
{
  ///
  /// determine number of particles from local rank to remote ranks
  ///
  _partsTo.resize(noRanks);
  _partsFrom.resize(noRanks);
  for (size_t i = 0; i < noDomains; i++)
    {
      _partsTo[ domainToRank[i] ] = _indices[i].size();
    }

  ///
  /// determine number of particles from remote ranks to local ranks
  ///
  barrier();
  MPI::COMM_WORLD.Alltoall(&_partsTo[0], 1, MPI::INT,
                           &_partsFrom[0], 1, MPI::INT);
  barrier();

  ///
  /// determine send and receive offsets. total offset is not added
  /// for send offset, as there is no total offset in the send buffer.
  ///
  _sendOffset[0] = 0;
  _recvOffset[0] = _totalOffset;
  for (size_t i = 1; i < noRanks; i++)
    {
      _sendOffset[i] = _sendOffset[i - 1] + _partsTo[i - 1];
      _recvOffset[i] = _recvOffset[i - 1] + _partsFrom[i - 1];
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

double CommunicationManager::wtime()
{
  return MPI::Wtime();
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
